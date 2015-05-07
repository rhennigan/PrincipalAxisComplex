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

GetPrincipalAxis     :: usage = ""
MakeTree             :: usage = ""
AverageClusterRadius :: usage = ""
ClusterDistance      :: usage = ""
ComputeTreeDistances :: usage = ""

Begin["`Private`"] (* Begin Private Context *)

$PathToJavaPlex = FileNameJoin[{NotebookDirectory[], "javaplex"}];
$PathToPerseus  = FileNameJoin[{NotebookDirectory[], "perseus"}];

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
  (* start direction can be chosen randomly, but this way it's deterministic *)
  (* this chooses the index of the coordinate with the most variance *)
  principalDimension = Last[Ordering[Variance[cluster]]];
  (* get the basis vector for the initial direction *)
  initialDirection = IdentityMatrix[dim][[principalDimension]];
  (* translates the cluster so that the centroid is on the origin *)
  shift = Shift[cluster];
  (* in practice, very few iterations are needed for a good split *)
  (* the 50 iteration limit is here to ensure completion if choosing random inital directions *)
  FixedPoint[Normalize[Total[shift.# * shift]] &, initialDirection, 50]
];

PartitionCluster[address_, cluster_List, limit_, p_] := Module[
  {len, shift, paxis, projection, order, data, separations,
    lowerBound, upperBound, largestGap, splitLocation,
    leftIndex, rightIndex,
    leftCluster, rightCluster,
    leftAddress, rightAddress},

(* memoization of cluster-address relationships *)
  GetCluster[address] = cluster;
  GetAddress[cluster] = address;

  len = Length[cluster];

  If[len < limit, (* maximum tree depth has been reached *)

    ((* then *)
      Sow[cluster]; (* keep this cluster as a leaf *)
      GoDownL[address] = address; (* going down when already at       *)
      GoDownR[address] = address; (*   the bottom doesn't go anywhere *)
    ),

    ((* else *)
      shift = Shift[cluster];
      paxis = GetPrincipalAxis[cluster];
      projection = shift.paxis; (* we're now working in a 1-dimensional space *)
      order = Ordering[projection];
      data = projection[[order]];

      (* find the largest gap to make a good hyperplane split *)
      separations = Differences[data];

      (* p is a balancing parameter:
         for p = 0, the split will be made wherever the gap is largest,
         which could result in a very unbalanced tree.
         for p = 1, the split always happens in the middle, so that both
         the left and right subtrees have an equal number of points*)
      lowerBound = Ceiling[p / 2 * len];
      upperBound = Ceiling[(1 - p / 2) * len];
      largestGap = Last[Ordering[separations]];
      splitLocation =
          If[lowerBound <= largestGap <= upperBound,
            largestGap,
            Ceiling[len / 2]
          ];
      leftIndex = order[[;; splitLocation]];
      rightIndex = order[[splitLocation + 1 ;;]];
      leftCluster = cluster[[leftIndex]];
      rightCluster = cluster[[rightIndex]];
      leftAddress = address * 2;
      rightAddress = address * 2 + 1;

      (* recursively descend into subtrees *)
      PartitionCluster[leftAddress, leftCluster, limit, p];
      PartitionCluster[rightAddress, rightCluster, limit, p];
    )
  ];

];

(* when a cluster has fewer than "sizeLimit" points, it will no longer be split *)
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

(* if no C compiler is available, change the "CompilationTarget -> "C"" option to "WVM" instead *)
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
            normalized = a / norm;
            axis = {#} & /@ normalized;
            axisT = Transpose[axis];
            projectedDown1 = shift1.axis;
            projectedDown2 = shift2.axis;
            p1 = Max[projectedDown1];
            p2 = Min[projectedDown2];

            errors1 =
                Table[Module[{projectionDistance, centerDistance, error},
                  projectionDistance = Norm[(v - offset) - (v - offset).axis.axisT];
                  centerDistance = Max[p1 - (v - offset).axis, (v - offset).axis - p2];
                  error = Norm[{centerDistance, projectionDistance}];
                  Prepend[v, error]], {v, trimmed1}];

            errors2 =
                Table[Module[{projectionDistance, centerDistance, error},
                  projectionDistance = Norm[(v - offset) - (v - offset).axis.axisT];
                  centerDistance = Max[p1 - (v - offset).axis, (v - offset).axis - p2];
                  error = Norm[{centerDistance, projectionDistance}];
                  Prepend[v, error]], {v, trimmed2}];

            len1 = Ceiling[Length[trimmed1] / 2];
            len2 = Ceiling[Length[trimmed2] / 2];
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

(* NeighborhoodScale controls the size of the neighborhood around a cluster, which is modified by the cluster's radius
   in order to adapt to local density changes *)
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
      keepQueue = Select[queue, ClusterDistance[current, #] < clusterDiameter * neighborhoodScale &];

      (* get left and right node numbers *)
      nextLeft = GoDownL[current];
      nextRight = GoDownR[current];

      (* check if left and right nodes exist *)
      leftExistsQ = nextLeft =!= current;
      rightExistsQ = nextRight =!= current;
      pairExistsQ = leftExistsQ && rightExistsQ;

      (* add them to the new queues if both exist *)
      leftQueue = If[pairExistsQ, Append[keepQueue, nextRight], keepQueue];
      rightQueue = If[pairExistsQ, Append[keepQueue, nextLeft], keepQueue];

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

(* computes the n-dimensional volume of the simplex defined by the n+1 vertices *)
SimplexContent[vertices_] := Module[
  {n, matrix, i, j},
  n = Length[vertices] - 1;
  matrix = ReplacePart[
    PadLeft[
      Table[
        Norm[vertices[[i]] - vertices[[j]]]^2,
        {i, n + 1}, {j, n + 1}
      ], {n + 2, n + 2}, 1
    ], {1, 1} -> 0
  ];
  Chop[Sqrt[((-1)^(n + 1) Det[matrix]) / (2^n (n!)^2)]]
]

RefineClustersGeneral[clusters_] := Module[
  {combined, center, shiftedClusters, pts, len},
  combined = Flatten[clusters, 1];
  center = Mean[combined];
  shiftedClusters = Map[# - center&, clusters, {2}];
  pts = Mean /@ shiftedClusters;
  len = Length[clusters];
  Table[
    Module[
      {pt1, pt2, line, cp, cpPoints, inner, axisDist, projError, error, split},
      pt1 = pts[[i]];
      pt2 = Mean[Delete[pts, i]];
      line = Normalize[pt2 - pt1];
      cp = shiftedClusters[[i]].line;
      cpPoints = # * line& /@ cp;
      inner = Max[cp];
      axisDist = Abs[cp - inner];
      projError = Norm /@ (shiftedClusters[[i]] - cpPoints);
      error = 2.0 axisDist + projError;
      split = Ceiling[Length[clusters[[i]]] / 2];
      clusters[[i, Ordering[error][[;; split]]]]
    ], {i, len}]
]

SpanningSimplexSequence[clusters_] := Module[
  {sequence},
  sequence = FixedPointList[RefineClustersGeneral, clusters];
  Map[Mean, sequence, {2}]
]



AverageClusterRadius[cluster_] := Module[
  {m},
  m = Mean[cluster];
  Mean[Norm[# - m]& /@ cluster]
]

GoUp[1] := 1
GoUp[address_] := Floor[address / 2]
GoDownL[address_] := 2address
GoDownR[address_] := 2address + 1
Transfer[address_] := GoDownL[GoUp[address]] /; OddQ[address]
Transfer[address_] := GoDownR[GoUp[address]] /; EvenQ[address]
CurrentLevel[address_] := Floor[Log[2, address] + 1]

NextUp[address_] := Transfer[GoUp[address]]

ParentNodes[address_] := FixedPointList[GoUp, address][[2 ;; -2]]

InitialSearchNodes[address_] := Module[
  {parentNodePath},
  parentNodePath = ParentNodes[address];
  Reverse[Complement[FixedPointList[NextUp, Transfer[address]], parentNodePath]]
]

VisitNode[cluster_, node_, maxDistance_] := Module[
  {targetCluster, distance},
  targetCluster = GetCluster[node];
  If[targetCluster != {},
    If[DistanceTable[cluster, targetCluster] == Infinity,
      distance = MinSpanningSimplexContent[{cluster, targetCluster}];
      DistanceTable[cluster, targetCluster] = distance;
      DistanceTable[targetCluster, cluster] = distance;,
      distance = DistanceTable[cluster, targetCluster];
    ];
    Sow[{GetAddress[cluster], node} -> distance];
    If[distance < maxDistance,
      VisitNode[cluster, GoDownL[node], maxDistance];
      VisitNode[cluster, GoDownR[node], maxDistance];
    ];
  ];
];

GetDistances[address_, p_] := Module[
  {cluster, rad},
  cluster = GetCluster[address];
  rad = p * AverageClusterRadius[cluster];
  First[Last[Reap[VisitNode[cluster, #, rad]& /@ InitialSearchNodes[address]]]]
]

GetDistances[address_, p_, globalRad_] := Module[
  {cluster, rad},
  cluster = GetCluster[address];
  rad = globalRad;
  First[Last[Reap[VisitNode[cluster, #, rad]& /@ InitialSearchNodes[address]]]]
]

TreeDistances[clusters_, p_] := Module[
  {addressList},
  addressList = GetAddress /@ clusters;
  Clear[DistanceTable];
  DistanceTable[_, _] = Infinity;
  Flatten[GetDistances[#, p]& /@ addressList, 1]
]

TreeDistances[clusters_, p_, globalRad_] := Module[
  {addressList},
  addressList = GetAddress /@ clusters;
  Clear[DistanceTable];
  DistanceTable[_, _] = Infinity;
  Flatten[GetDistances[#, p, globalRad]& /@ addressList, 1]
]

DistanceMapper[cluster_, rad_, target_] := Module[
  {distance},
  distance = MinSpanningSimplexContent[{cluster, target}];
  If[distance < rad, {distance, {GoDownL[node], GoDownR[node]}}, {distance, {}}]
]

TreeSearch[address_, queue_, leaves_, globalRad_] := Module[
  {cluster, distances, filteredDistances, carry},
  cluster = GetCluster[address];
  If[cluster != {} && queue != {},
    distances = {#, MinSpanningSimplexContent[{cluster, GetCluster[#]}]}& /@ queue;
    filteredDistances = Select[distances, #[[2]] < globalRad&];
    Sow[{address, #[[1]]} -> #[[2]]]& /@ filteredDistances;
    carry = Union[Flatten[{GoDownL[#], GoDownR[#]}& /@ filteredDistances[[All, 1]]]];
    If[Not[MemberQ[leaves, address]],
      TreeSearch[GoDownL[address], Prepend[carry, GoDownR[address]], leaves, globalRad];
      TreeSearch[GoDownR[address], Prepend[carry, GoDownL[address]], leaves, globalRad];
    ];
  ]
]

TreeSearch[clusters_, globalRad_] := Module[
  {leaves},
  leaves = GetAddress /@ clusters;
  First[Last[Reap[TreeSearch[2, {3}, leaves, globalRad]]]]
]

TreeSearch[address_, queue_, p_] := Module[
  {rad, distances, filter, nextQueue, left, right},

  rad = p * AverageClusterRadius[GetCluster[address]];
  distances = {#, MinSpanningSimplexContent[address, #]}& /@ queue;

  (*len=Length[GetCluster[address]];*)

  Scan[(
    Sow[{{address, #[[1]]}, {Length[GetCluster[address]], Length[GetCluster[#[[1]]]]}}];
    DistanceTable[address, #[[1]]] = #[[2]];
    DistanceTable[#[[1]], address] = #[[2]];
  )&, distances];
  filter = Select[distances, #[[2]] < rad&];

  nextQueue = Union[Flatten[{GoDownL[#], GoDownR[#]}& /@ Prepend[filter[[All, 1]], address], 1]];

  left = GoDownL[address];
  right = GoDownR[address];

  (*test={left,right,address,nextQueue,queue};
Print[{"Left: "<>ToString[left]<>"   Right: "<>ToString[right]<>"   Address: "<>ToString[address]<>"   Next: "<>ToString[nextQueue]<>"   Queue: "<>ToString[queue]}];*)
  If[queue != {} || address == 1,

    If[Not[left == address && nextQueue == queue],
      TreeSearch[left, Complement[nextQueue, queue], p];
    ];

    If[Not[right == address && nextQueue == queue],
      TreeSearch[right, Complement[nextQueue, queue], p];
    ];
  ];
]

Clear[MinSpanningSimplexContent]

MinSpanningSimplexContent[clusters_] := Module[
  {d, scale},
  If[Not[MemberQ[clusters, {}]],
    d = Length[clusters];
    scale = Sqrt[2^(-1 + d) / d] Gamma[d];
    scale * Min[SimplexContent /@ SpanningSimplexSequence[clusters]],
    Infinity]
]


MinSpanningSimplexContent[cluster1_, cluster2_] := MinSpanningSimplexContent[{cluster1, cluster2}]
MinSpanningSimplexContent[a_Integer, b_Integer] := (
  If[
    DistanceTable[a, b] == Infinity,
    DistanceTable[a, b] = MinSpanningSimplexContent[GetCluster[a], GetCluster[b]];
    DistanceTable[b, a] = DistanceTable[a, b];
  ];DistanceTable[a, b]) /; a != b
MinSpanningSimplexContent[a_Integer, a_Integer] := 0

FilterCluster[cluster_, p_] := Module[
  {order},
  order = Ordering[Norm /@ Shift[cluster]];
  cluster[[order]][[;; Ceiling[p * Length[cluster]]]]
]

ParallelTreeSearch[address_, queue_, p_] := Module[
  {rad, distances, filter, nextQueue, left, right},

  rad = p * AverageClusterRadius[GetCluster[address]];
  distances = {#, ParallelContent[address, #]}& /@ queue;

  filter = Select[distances, #[[2]] < rad&];

  Sow[{{address, #[[1]]}, #[[2]]}]& /@ filter;
  Sow[{{#[[1]], address}, #[[2]]}]& /@ filter;

  nextQueue = Union[Flatten[{GoDownL[#], GoDownR[#]}& /@ Prepend[filter[[All, 1]], address], 1]];

  left = GoDownL[address];
  right = GoDownR[address];

  If[queue != {} || address == 1,

    If[Not[left == address && nextQueue == queue],
      ParallelTreeSearch[left, Complement[nextQueue, queue], p];
    ];

    If[Not[right == address && nextQueue == queue],
      ParallelTreeSearch[right, Complement[nextQueue, queue], p];
    ];
  ];
]

ParallelTreeSearch[address_, queue_, p_, q_] := Module[
  {rad, distances, filter, nextQueue, left, right},

  rad = p * AverageClusterRadius[GetCluster[address]] / q;
  distances = {#, ParallelContent[address, #, q]}& /@ queue;

  filter = Select[distances, #[[2]] < rad&];

  Sow[{{address, #[[1]]}, #[[2]]}]& /@ filter;
  Sow[{{#[[1]], address}, #[[2]]}]& /@ filter;

  nextQueue = Union[Flatten[{GoDownL[#], GoDownR[#]}& /@ Prepend[filter[[All, 1]], address], 1]];

  left = GoDownL[address];
  right = GoDownR[address];

  If[queue != {} || address == 1,

    If[Not[left == address && nextQueue == queue],
      ParallelTreeSearch[left, Complement[nextQueue, queue], p, q];
    ];

    If[Not[right == address && nextQueue == queue],
      ParallelTreeSearch[right, Complement[nextQueue, queue], p, q];
    ];
  ];
]

ParallelContent[cluster1_, cluster2_] := MinSpanningSimplexContent[{cluster1, cluster2}]
ParallelContent[a_Integer, b_Integer] := MinSpanningSimplexContent[GetCluster[a], GetCluster[b]] /; a != b
ParallelContent[a_Integer, a_Integer] := 0

ParallelContent[cluster1_, cluster2_, p_] := MinSpanningSimplexContent[{FilterCluster[cluster1, p], FilterCluster[cluster2, p]}]
ParallelContent[a_Integer, b_Integer, p_] := MinSpanningSimplexContent[FilterCluster[GetCluster[a], p], FilterCluster[GetCluster[b], p]] /; a != b
ParallelContent[a_Integer, a_Integer, p_] := 0

CustomColor[i_] := {
  RGBColor[4 / 85, 6 / 17, 43 / 85],
  RGBColor[214 / 255, 74 / 255, 13 / 85],
  RGBColor[2 / 5, 203 / 255, 19 / 85],
  RGBColor[142 / 255, 173 / 255, 57 / 85],
  RGBColor[3 / 5, 52 / 255, 66 / 85],
  RGBColor[43 / 51, 224 / 255, 229 / 255]
}[[i]]

\[CapitalGamma][x_] := Gamma[x]
SphereVolume[r_, d_] := (\[Pi]^(d / 2) r^d) / \[CapitalGamma][1 + d / 2]
ClusterDensity[cluster_] := 0.0 /; Length[cluster] == 1
ClusterDensity[cluster_] := Length[cluster] / SphereVolume[AverageClusterRadius[cluster], Length[First[cluster]]]

CutOutliers[clusters_, p_, size_] := Module[
  {densities, q, lowDensity},
  densities = ClusterDensity /@ clusters;
  q = Quantile[densities, p];
  lowDensity = Select[clusters, ClusterDensity[#] < q || Length[#] <= size&];
  Sow[GetAddress[#]]& /@ lowDensity;
  Scan[Function[cluster, Module[
    {address, parent},
    address = GetAddress[cluster];
    parent = GoUp[address];
    GoDownL[parent] = parent;
    GoDownR[parent] = parent;
  ]], lowDensity]
]

GetTrianglesMaybe[i_, g_, dim_] := With[{ng = NeighborhoodGraph[g, i]}, Sort /@ FindClique[ng, dim + 1, All]]

GetTrianglesMaybe[i_, g_] := With[{ng = NeighborhoodGraph[g, i]}, Sort /@ FindClique[ng, Infinity, All]]

RecordSimplexAge[simplex_, clusters_] := If[SimplexAge[simplex] == -1,
  Module[
    {lowerLimit},
    lowerLimit = Max[SimplexAge /@ Subsets[simplex, {Length[simplex] - 1}]];
    If[lowerLimit == Infinity,
      SimplexAge[simplex] = Infinity,
      SimplexAge[simplex] = lowerLimit + MinSpanningSimplexContent[clusters[[simplex]]]
    ];
    Sow[{simplex, SimplexAge[simplex]}];
  ]
]

RecordSimplexAgeVR[simplex_, clusters_] := If[SimplexAge[simplex] == -1,
  SimplexAge[simplex] = Max[SimplexAge /@ Subsets[simplex, {Length[simplex] - 1}]];
  Sow[{simplex, SimplexAge[simplex]}];
]

AgeFromIndex[index_] := If[index == Infinity, Infinity, filteredComplex[[index, 2]]]

Needs["JLink`"];
With[
  {pathToJavaPlex = $PathToJavaPlex},
  AddToClassPath[pathToJavaPlex];
  InstallJava[];
  ReinstallJava[JVMArguments -> "-Xmx16384m"];
  LoadJavaClass["edu.stanford.math.plex4.api.Plex4"];
  apiPlex4 = JavaNew["edu.stanford.math.plex4.api.Plex4"];
  LoadJavaClass["edu.stanford.math.plex4.examples.PointCloudExamples"];
  examplesPointCloudExamples = JavaNew["edu.stanford.math.plex4.examples.PointCloudExamples"];
  LoadJavaClass["edu.stanford.math.plex4.metric.impl.ExplicitMetricSpace"];
(*ExplicitMetricSpace=JavaNew["edu.stanford.math.plex4.metricimpl.ExplicitMetricSpace"];*)
];

End[] (* End Private Context *)

EndPackage[]