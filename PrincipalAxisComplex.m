(* Mathematica Package                   *)
(* Created by IntelliJ IDEA              *)

(* :Title: PrincipalAxisComplex          *)
(* :Context: PrincipalAxisComplex`       *)
(* :Author: Richard Hennigan             *)
(* :Date: 5/1/2015                       *)

(* :Package Version: 1.0                 *)
(* :Mathematica Version:                 *)
(* :Copyright: (c) 2015 Richard Hennigan *)
(* :Keywords:                            *)
(* :Discussion:                          *)

BeginPackage["PrincipalAxisComplex`"]
(* Exported symbols added here with SymbolName::usage *)

Begin["`Private`"] (* Begin Private Context *)

Reset[] := (
  Clear[GetCluster, GetAddress, GetPoint];

  (* default values *)
  GetCluster[___] := {};
  GetAddress[___] := -1;
  GetPoint[___] := {};
  GoDownL[address_] := 2 address;
  GoDownR[address_] := 2 address + 1;
)

Shift[cluster_] := Module[{m = Mean[cluster]}, # - m & /@ cluster];

GetPrincipalAxis[cluster_] := Module[
  {dim, principalDimension, shift, initialDirection},
  dim = Length[First[cluster]];
  principalDimension = Last[Ordering[Variance[cluster]]];
  shift = Shift[cluster];
  initialDirection = IdentityMatrix[dim][[principalDimension]];
  FixedPoint[Normalize[Total[shift.#*shift]] &, initialDirection, 50]
];

PartitionCluster[address_, cluster_List, limit_, p_] := Module[
  {len, shift, paxis, projection, order, data, separations,
    lowerBound, upperBound, largestGap, splitLocation,
    leftIndex, rightIndex,
    leftCluster, rightCluster,
    leftAddress, rightAddress},

  GetCluster[address] = cluster;
  GetAddress[cluster] = address;

  len = Length[cluster];

  If[len < limit,

    Sow[cluster];
    GoDownL[address] = address;
    GoDownR[address] = address;

    ,

    shift = Shift[cluster];
    paxis = GetPrincipalAxis[cluster];
    projection = shift.paxis;
    order = Ordering[projection];
    data = projection[[order]];
    separations = Differences[data];

    lowerBound = Ceiling[p/2*len];
    upperBound = Ceiling[(1 - p/2)*len];
    largestGap = Last[Ordering[separations]];
    splitLocation =
        If[lowerBound <= largestGap <= upperBound, largestGap,
          Ceiling[len/2]
        ];

    leftIndex = order[[;; splitLocation]];
    rightIndex = order[[splitLocation + 1 ;;]];

    leftCluster = cluster[[leftIndex]];
    rightCluster = cluster[[rightIndex]];

    leftAddress = address*2;
    rightAddress = address*2 + 1;

    PartitionCluster[leftAddress, leftCluster, limit, p];
    PartitionCluster[rightAddress, rightCluster, limit, p];

  ];

];

MakeTree[points_, sizeLimit_, p_] := (
  Clear[GetCluster, GetAddress];
  GetCluster[_] = {};
  GetAddress[_] = -1;
  First[Last[Reap[PartitionCluster[1, points, sizeLimit, p]]]]
);

AverageClusterRadius[cluster_] := Module[
  {m},
  m = Mean[cluster];
  Mean[Norm[# - m] & /@ cluster]
];

ClusterDistanceC =
    Compile[{{cluster1, _Real, 2}, {cluster2, _Real, 2}},
      If[Norm[Mean[cluster1] - Mean[cluster2]] < 0.000001, 0.0,
        Module[
          {dim, offset, shift1, shift2, m1, m2, a, norm, normalized, axis,
            axisT, projectedDown1, projectedDown2, p1, p2, errors1, errors2,
            len1, len2, trimmed1, trimmed2},

          dim = Dimensions[cluster1][[2]];

          trimmed1 = cluster1;
          trimmed2 = cluster2;

          While[Length[trimmed1] > 1 || Length[trimmed2] > 1,
            offset = Mean[Join[trimmed1, trimmed2]];
            shift1 = # - offset & /@ trimmed1;
            shift2 = # - offset & /@ trimmed2;
            m1 = Mean[shift1];
            m2 = Mean[shift2];
            a = m2 - m1;
            norm = Norm[a];
            normalized = a/norm;
            axis = {#} & /@ normalized;
            axisT = Transpose[axis];
            projectedDown1 = shift1.axis;
            projectedDown2 = shift2.axis;
            p1 = Max[projectedDown1];
            p2 = Min[projectedDown2];

            errors1 =
                Table[Module[{projectionDistance, centerDistance, error},
                  projectionDistance =
                      Norm[(v - offset) - (v - offset).axis.axisT];
                  centerDistance =
                      Max[p1 - (v - offset).axis, (v - offset).axis - p2];
                  error = Norm[{centerDistance, projectionDistance}];
                  Prepend[v, error]
                ], {v, trimmed1}];

            errors2 =
                Table[Module[{projectionDistance, centerDistance, error},
                  projectionDistance =
                      Norm[(v - offset) - (v - offset).axis.axisT];
                  centerDistance =
                      Max[p1 - (v - offset).axis, (v - offset).axis - p2];
                  error = Norm[{centerDistance, projectionDistance}];
                  Prepend[v, error]
                ], {v, trimmed2}];

            len1 = Ceiling[Length[trimmed1]/2];
            len2 = Ceiling[Length[trimmed2]/2];
            trimmed1 = SortBy[errors1, First][[;; len1, 2 ;; dim + 1]];
            trimmed2 = SortBy[errors2, First][[;; len2, 2 ;; dim + 1]];
          ];

          Norm[First[trimmed1] - First[trimmed2]]
        ]
      ],
      CompilationTarget -> "C",
      RuntimeAttributes -> {Listable},
      Parallelization -> True,
      RuntimeOptions -> "Speed"
    ];

Options[ComputeTreeDistances] = {"NeighborhoodScale" -> 0.5};
ComputeTreeDistances[current_, queue_, opts : OptionsPattern[]] :=
    Module[
      {
        cluster, clusterDiameter, neighborhoodScale,
        queuePaired, keepQueue,
        nextLeft, nextRight,
        leftExistsQ, rightExistsQ, pairExistsQ,
        leftQueue, rightQueue
      },

      neighborhoodScale = OptionValue["NeighborhoodScale"];
      cluster = GetCluster[current];
      clusterDiameter = 2 AverageClusterRadius[cluster];
      queuePaired = {current, #} & /@ queue;

      (* compute and store distances *)
      ClusterDistance @@@ queuePaired;
      Sow /@ queuePaired;

      (* keep close neighbors *)
      keepQueue =
          Select[queue,
            ClusterDistance[current, #] <
                clusterDiameter*neighborhoodScale &];

      (* get left and right node numbers *)
      nextLeft = GoDownL[current];
      nextRight = GoDownR[current];

      (* check if left and right nodes exist *)
      leftExistsQ = nextLeft =!= current;
      rightExistsQ = nextRight =!= current;
      pairExistsQ = leftExistsQ && rightExistsQ;

      (* add them to the new queues if both exist *)
      leftQueue =
          If[pairExistsQ, Append[keepQueue, nextRight], keepQueue];
      rightQueue =
          If[pairExistsQ, Append[keepQueue, nextLeft], keepQueue];

      (* recursively descend if child nodes exist *)
      If[leftExistsQ, ComputeTreeDistances[nextLeft, leftQueue, opts]];
      If[rightExistsQ, ComputeTreeDistances[nextRight, rightQueue, opts]];
    ];

ClearAll[ClusterDistance];
ClusterDistance[___] := Infinity;
ClusterDistance[q1_Integer, q2_Integer] /; q1 === q2 := 0.0;
ClusterDistance[q1_Integer, q2_Integer] := Module[
  {distance = ClusterDistanceC[GetCluster[q1], GetCluster[q2]]},
  ClusterDistance[q1, q2] = ClusterDistance[q2, q1] = distance;
];

End[] (* End Private Context *)

EndPackage[]