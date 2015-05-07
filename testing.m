(* Mathematica source file  *)
(* Created by IntelliJ IDEA *)
(* :Author: rhennigan *)
(* :Date: 5/4/2015 *)

Needs["PrincipalAxisComplex`"]
Needs["HypothesisTesting`"]

selectGalaxy[img_, t_] := SelectComponents[FillingTransform @ Binarize[img, t], "Area", -1]
selectGalaxy[img_] := SelectComponents[FillingTransform @ Binarize[img], "Area", -1]
galaxyAngle[galaxyComponent_] := 1 /. ComponentMeasurements[galaxyComponent, "Orientation"]
croppedCentroid[img_, s_] := ImageMeasurements[ImageCrop[img, s], "IntensityCentroid"] / s
contrast[img_, exponent_] := ImageApply[Max[0.0, #]^exponent &, img]

selectCenterGalaxy[img_, t_] :=
    Module[
      {components, label},
      components = FillingTransform @ Binarize[img, t];
      label = MinimalBy[ComponentMeasurements[components, "Centroid"], Norm[#[[2]] - ImageDimensions[img] / 2] &][[1, 1]];
      SelectComponents[components, "Label", # == label &]
    ]

selectCenterGalaxy[img_] :=
    Module[
      {components, label},
      components = FillingTransform @ Binarize[img];
      label = MinimalBy[ComponentMeasurements[components, "Centroid"], Norm[#[[2]] - ImageDimensions[img] / 2] &][[1, 1]];
      SelectComponents[components, "Label", # == label &]
    ]

transformImage[img_] :=
    Module[{w, h, components, label, img2, imgc, orientation, rotated,
      cropped, resized},
      {w, h} = ImageDimensions[img];
      components = DeleteSmallComponents[FillingTransform @ Dilation[AlphaChannel[RemoveBackground[img, {{{1, 1}, {w, 1}, {1, h}, {w, h}}, .025}]], 2]];
      label = MinimalBy[ComponentMeasurements[components, "Centroid"], Norm[#[[2]] - ImageDimensions[img] / 2] &][[1, 1]];
      img2 = ImageMultiply[img, Blur[SelectComponents[components, "Label", # == label &], 8]];
      imgc = ImageCrop[img2, {w, h}];
      orientation = galaxyAngle[selectGalaxy[imgc]];
      rotated = ImageCrop[ImageRotate[imgc, -orientation], Min[w, h] / 2];
      cropped = ImageTrim[rotated, {{-20, -20}, {20, 20}} + Round[1 /. ComponentMeasurements[selectCenterGalaxy[rotated, .05], "BoundingBox"]]];
      resized = With[{m = Max[ImageDimensions[cropped]]}, ImageResize[ImageCrop[cropped, {m, m}], Min[w, h] / 2]];
      resized
    ]

processedImageDir = "E:\\Mathematica Frames\\filteredGalaxies";
originalImageDir = "E:\\GDrive\\Classes\\Machine Learning\\project\\data\\images_training_rev1\\images_training_rev1";

With[{import = Import["E:\\GDrive\\Classes\\Machine Learning\\project\\data\\training_solutions_rev1\\training_solutions_rev1.csv"]},
  variableNames = First[import];
  data = Rest[import];
];
data[[All, 2 ;;]] = N[data[[All, 2 ;;]]];

idHasProcessedFile[id_] := FileExistsQ[FileNameJoin[{processedImageDir, ToString[id] <> "_out.jpg"}]]
entryHasProcessedFile[entry_] := idHasProcessedFile[First[entry]]
(*getProcessedFile[id_]:=Import[FileNameJoin[{processedImageDir, ToString[id]<>"_out.jpg"}]]*)

loadStruct[entry_] :=
    Module[{id, assoc},
      id = First[entry];
      assoc = <|
          "smooth" -> <|
              "smooth" -> entry[[2]],
              "featuresOrDisk" -> entry[[3]],
              "starOrArtifact" -> entry[[4]]
              |>,
          "edgeOn" -> <|
              "true" -> entry[[5]],
              "false" -> entry[[6]]
              |>,
          "centerBar" -> <|
              "true" -> entry[[7]],
              "false" -> entry[[8]]
              |>,
          "spiralArm" -> <|
              "true" -> entry[[9]],
              "false" -> entry[[10]]
              |>,
          "bulge" -> <|
              "a1" -> entry[[11]],
              "a2" -> entry[[12]],
              "a3" -> entry[[13]],
              "a4" -> entry[[14]]
              |>,
          "odd" -> <|
              "a1" -> entry[[15]],
              "a2" -> entry[[16]]
              |>,
          "round" -> <|
              "a1" -> entry[[17]],
              "a2" -> entry[[18]],
              "a3" -> entry[[19]]
              |>,
          "oddFeature" -> <|
              "a1" -> entry[[20]],
              "a2" -> entry[[21]],
              "a3" -> entry[[22]],
              "a4" -> entry[[23]],
              "a5" -> entry[[24]],
              "a6" -> entry[[25]],
              "a7" -> entry[[26]]
              |>,
          "bulgeShape" -> <|
              "a1" -> entry[[27]],
              "a2" -> entry[[28]],
              "a3" -> entry[[29]]
              |>,
          "tightlyWoundArms" -> <|
              "a1" -> entry[[30]],
              "a2" -> entry[[31]],
              "a3" -> entry[[32]]
              |>,
          "numberOfArms" -> <|
              "a1" -> entry[[33]],
              "a2" -> entry[[34]],
              "a3" -> entry[[35]],
              "a4" -> entry[[36]],
              "a5" -> entry[[37]],
              "a6" -> entry[[38]]
              |>
          |>;
      id -> assoc
    ]
sigmoid[x_] := 1 / (1 + E^-x)

inverseSigmoid[y0_] := With[{y = $MachineEpsilon + y0 - 2 $MachineEpsilon y0}, Log[-(y / (-1 + y))]]

taskRelationSmooth[featureFunction_, id_] :=
    Module[{img, features, smooth, featuresOrDisk, starOrArtifact},
      img = getProcessedFile[id];
      features = featureFunction[img];
      smooth = inverseSigmoid[lookup[id]["smooth"]["smooth"]];
      featuresOrDisk = inverseSigmoid[lookup[id]["smooth"]["featuresOrDisk"]];
      starOrArtifact = inverseSigmoid[lookup[id]["smooth"]["starOrArtifact"]];
      features -> {smooth, featuresOrDisk, starOrArtifact}
    ]

lookup = Association[loadStruct /@ data];
(*getProcessedFile=Association[ParallelMap[#\[Rule]Import[FileNameJoin[{processedImageDir,ToString[#]<>"_out.jpg"}]]&,Keys[lookup]]];*)
DumpSave["E:\\GDrive\\Classes\\Machine Learning\\project\\data\\lookup.mx", lookup];
DumpSave["E:\\GDrive\\Classes\\Machine Learning\\project\\data\\getProcessedFile.mx", getProcessedFile];

trueClassification[id_, task_] := Last[SortBy[Keys[lookup[id][task]], lookup[id][task]]]

testClassifier[featureFunction_, c_, id_] :=
    Module[{img},
      img = featureFunction[getProcessedFile[id]];
      {c[img], trueClassification[id, "smooth"]}
    ]

assesFeatureFunction[featureFunction_, trainingSetSize_, sampleCount_, sampleSize_, method_, parallel_] :=
    Module[
      {
        time, dataSample, lookupSample, linked, c,
        validationSampleKeys, validationSampleValues, validationSampleImages, validationSamples,
        correct, ci, mean
      },
      time = First[AbsoluteTiming[
        dataSample = RandomSample[data, trainingSetSize];
        lookupSample = Association[loadStruct /@ dataSample];

        linked = featureFunction[getProcessedFile[#]] ->
            Last[SortBy[Keys[lookupSample[#]["smooth"]], lookupSample[#]["smooth"]]]& /@ Keys[lookupSample];

        c = Classify[linked, Method -> method];

        validationSampleKeys = Table[RandomSample[Keys[lookup], sampleSize], {sampleCount}];
        validationSampleValues = Map[trueClassification[#, "smooth"]&, validationSampleKeys, {2}];
        validationSampleImages = Map[getProcessedFile, validationSampleKeys, {2}];
        validationSamples = Transpose[{validationSampleKeys, validationSampleValues, validationSampleImages}, {3, 1, 2}];

        correct = Flatten[If[parallel, ParallelTable, Table][Module[
          {validationSampleIDs, predictedValues, trueValues, images, correctCount},
          validationSampleIDs = validationSample[[All, 1]];
          trueValues = validationSample[[All, 2]];
          images = validationSample[[All, 3]];
          predictedValues = c[featureFunction[#]]& /@ images;
          correctCount = Boole[SameQ @@@ Transpose[{predictedValues, trueValues}]]
        (*N[correctCount/Length[validationSampleIDs]]*)
        ], {validationSample, validationSamples}]];

        ci = MeanCI[correct, ConfidenceLevel -> 0.99];
        mean = Mean[N[correct]];
      ]];
      {c, mean, ci, time}
    ]

(* feature functions *)
luminance[img_] := ImageAdjust[ColorSeparate[ColorConvert[img, "LUV"]][[1]]]
chrominance1[img_] := ImageAdjust[ColorSeparate[ColorConvert[img, "LUV"]][[2]]]
chrominance2[img_] := ImageAdjust[ColorSeparate[ColorConvert[img, "LUV"]][[3]]]
grayscale[img_] := ColorConvert[img, "Grayscale"]
gradient1[img_] := GradientFilter[img, 1]
gradient2[img_] := GradientFilter[img, 2]
gradient3[img_] := GradientFilter[img, 3]
gradient4[img_] := GradientFilter[img, 4]
red[img_] := ImageAdjust[ColorSeparate[img][[1]]]
green[img_] := ImageAdjust[ColorSeparate[img][[2]]]
blue[img_] := ImageAdjust[ColorSeparate[img][[3]]]

kldGradients[img_] := ColorCombine[(ImageAdjust /@ KarhunenLoeveDecomposition[ColorConvert[#, "Grayscale"]& /@ Join[ColorSeparate[ColorConvert[img, "LUV"]], Table[ImageAdjust[GradientFilter[img, i]], {i, 1, 8, 3}]]][[1]])[[;; 3]]]

functions = {luminance, chrominance1, chrominance2, grayscale, red, green, blue, gradient1, gradient2, gradient3, gradient4};
methods = {"LogisticRegression", "NaiveBayes", "NearestNeighbors", "RandomForest", "SupportVectorMachine"};

trainingSetSize = 200;
sampleCount = 8; (* parallelized parameter *)
sampleSize = 100;

testResults = Table[{m, ToString[f], assesFeatureFunction[f, trainingSetSize, sampleCount, sampleSize, m, True]}, {m, methods}, {f, functions}];
testResultsF = Flatten[testResults, 1];

Row[{
(* method accuracy *)
  TableForm[Reverse[SortBy[{#, Function[m, Mean[Select[testResultsF, #[[1]] == ToString[m]&][[All, 3, 2]]]][#]}& /@ methods, Last]], TableHeadings -> {None, {"Method", "Accuracy"}}],
  Spacer[10],
(* method time *)
  TableForm[SortBy[{#, Function[m, Mean[Select[testResultsF, #[[1]] == ToString[m]&][[All, 3, 4]]]][#]}& /@ methods, Last], TableHeadings -> {None, {"Method", "Time (s)"}}]
}]

Row[{
(* function accuracy *)
  TableForm[Reverse[SortBy[{#, Function[f, Mean[Select[testResultsF, #[[2]] == ToString[f]&][[All, 3, 2]]]][#]}& /@ functions, Last]], TableHeadings -> {None, {"Function", "Accuracy"}}],
  Spacer[10],
(* function time *)
  TableForm[SortBy[{#, Function[f, Mean[Select[testResultsF, #[[2]] == ToString[f]&][[All, 3, 4]]]][#]}& /@ functions, Last], TableHeadings -> {None, {"Function", "Time (s)"}}]
}]

persistence[steps_, featureFunction_, img_] := Module[
  {k, w, h, cwComplexData, cwComplexString, imgOut0, imgOut1, barcodes0, subsetBarcodes0, barcodes1, subsetBarcodes1},
  k = ToString[$KernelID];
  {w, h} = ImageDimensions[img];
  cwComplexData = Flatten[Floor[1 + steps - steps * Rescale[ImageData[ColorConvert[featureFunction[img], "Grayscale"]]]]];
  cwComplexString = StringJoin[
    "2", "\n",
    ToString[w], "\n",
    ToString[h], "\n",
    Riffle[ToString /@ cwComplexData, "\n"]
  ];
  Export[FileNameJoin[{$PathToPerseus, "perseusInput" <> k}], cwComplexString, "TEXT"];
  SetDirectory[$PathToPerseus];
  Get["!perseusWin.exe cubtop perseusInput" <> k <> " perseusOutput" <> k] // Quiet;
  imgOut0 = ConstantArray[0, {h, w}];
  imgOut1 = ConstantArray[0, {h, w}];
  barcodes0 = Partition[ToExpression /@ StringSplit[Import[FileNameJoin[{$PathToPerseus, "perseusOutput" <> k <> "_0.txt"}]]], 2] /. (-1 -> steps);
  subsetBarcodes0 = SortBy[barcodes0, #[[1]] - #[[2]]&][[;; Min[steps, Length[barcodes0]]]];
  barcodes1 = Partition[ToExpression /@ StringSplit[Import[FileNameJoin[{$PathToPerseus, "perseusOutput" <> k <> "_1.txt"}]]], 2] /. (-1 -> steps);
  subsetBarcodes1 = SortBy[barcodes1, #[[1]] - #[[2]]&][[;; Min[steps, Length[barcodes1]]]];

  Do[Module[{b1, b2, len},
    {b1, b2} = subsetBarcodes0[[i]];
    len = b2 - b1;
    imgOut0[[i, Max[1, b1] ;; Min[w, b2]]] = 1
  ], {i, Length[subsetBarcodes0]}];

  Do[Module[{b1, b2, len},
    {b1, b2} = subsetBarcodes1[[i]];
    len = b2 - b1;
    imgOut1[[Min[h, Max[1, h - i + 1]], Max[1, w - b2 + 1] ;; Min[w, w - b1 + 1]]] = 1
  ], {i, Length[subsetBarcodes1]}];

  ColorCombine[{Image[imgOut0], Image[imgOut1]}]
]

persistence2[steps_, featureFunction_, img_] := Module[
  {k, w, h, src, cwComplexData, cwComplexString, barcodes0, barcodes1, arr, i0, i1},
  k = ToString[$KernelID];
  {w, h} = ImageDimensions[img];
  src = ColorConvert[featureFunction[img], "Grayscale"];
  cwComplexData = Flatten[Floor[1 + steps - steps * Rescale[ImageData[src]]]];
  cwComplexString = StringJoin[
    "2", "\n",
    ToString[w], "\n",
    ToString[h], "\n",
    Riffle[ToString /@ cwComplexData, "\n"]
  ];
  Export[FileNameJoin[{$PathToPerseus, "perseusInput" <> k}], cwComplexString, "TEXT"];
  SetDirectory[$PathToPerseus];
  Get["!perseusWin.exe cubtop perseusInput" <> k <> " perseusOutput" <> k] // Quiet;
  barcodes0 = Partition[ToExpression /@ StringSplit[Import[FileNameJoin[{$PathToPerseus, "perseusOutput" <> k <> "_0.txt"}]]], 2] /. (-1 -> w);
  barcodes1 = Partition[ToExpression /@ StringSplit[Import[FileNameJoin[{$PathToPerseus, "perseusOutput" <> k <> "_1.txt"}]]], 2] /. (-1 -> w);
  arr = ConstantArray[0, {64, 64}];
  i0 = N[Fold[Function[{a, b}, ReplacePart[a, {i_, j_} /; b[[1]] <= i <= b[[2]] && b[[1]] <= j <= b[[2]] :> (a[[i, j]] + b[[2]] - b[[1]])]], arr, barcodes0]];
  i1 = N[Fold[Function[{a, b}, ReplacePart[a, {i_, j_} /; b[[1]] <= i <= b[[2]] && b[[1]] <= j <= b[[2]] :> (a[[i, j]] + b[[2]] - b[[1]])]], arr, barcodes1]];
  ColorCombine[{Image[i0 / Max[i0, i1]], Image[i1 / Max[i0, i1]], src}]
];

(* feature functions *)
identity = Identity;
luma[img_] := ImageAdjust[ColorSeparate[ColorConvert[img, "LUV"]][[1]]]
chrominance1[img_] := ImageAdjust[ColorSeparate[ColorConvert[img, "LUV"]][[2]]]
chrominance2[img_] := ImageAdjust[ColorSeparate[ColorConvert[img, "LUV"]][[3]]]
grayscale[img_] := ColorConvert[img, "Grayscale"]
gradient1[img_] := ImageAdjust[GradientFilter[img, 1]]
gradient2[img_] := ImageAdjust[GradientFilter[img, 2]]
gradient3[img_] := ImageAdjust[GradientFilter[img, 3]]
gradient4[img_] := ImageAdjust[GradientFilter[img, 4]]
red[img_] := ImageAdjust[ColorSeparate[img][[1]]]
green[img_] := ImageAdjust[ColorSeparate[img][[2]]]
blue[img_] := ImageAdjust[ColorSeparate[img][[3]]]
gabor2x[img_] := ImageAdjust[GaborFilter[img, 2, {1, 0}]]
gabor2y[img_] := ImageAdjust[GaborFilter[img, 2, {0, 1}]]
gabor3x[img_] := ImageAdjust[GaborFilter[img, 3, {1, 0}]]
gabor3y[img_] := ImageAdjust[GaborFilter[img, 3, {0, 1}]]
gabor4x[img_] := ImageAdjust[GaborFilter[img, 4, {1, 0}]]
gabor4y[img_] := ImageAdjust[GaborFilter[img, 4, {0, 1}]]
gabor8x[img_] := ImageAdjust[GaborFilter[img, 8, {1, 0}]]
gabor8y[img_] := ImageAdjust[GaborFilter[img, 8, {0, 1}]]
corner2[img_] := ImageAdjust[CornerFilter[img, 2]]
corner3[img_] := ImageAdjust[CornerFilter[img, 3]]
corner4[img_] := ImageAdjust[CornerFilter[img, 4]]
ridge2[img_] := ImageAdjust[RidgeFilter[img, 2]]
ridge3[img_] := ImageAdjust[RidgeFilter[img, 3]]
ridge4[img_] := ImageAdjust[RidgeFilter[img, 4]]
laplacianGaussianFilter1[img_] := ImageAdjust[LaplacianGaussianFilter[img, 1]]
laplacianGaussianFilter2[img_] := ImageAdjust[LaplacianGaussianFilter[img, 2]]
laplacianGaussianFilter3[img_] := ImageAdjust[LaplacianGaussianFilter[img, 3]]
laplacian1[img_] := ImageAdjust[LaplacianFilter[img, 1]]
laplacian2[img_] := ImageAdjust[LaplacianFilter[img, 2]]

binList10[img_] := Module[{gimg = gradient3[img]},
  Table[DeleteSmallComponents[Erosion[Binarize[gimg, FindThreshold[img, Method -> {"BlackFraction", t}]], 1]], {t, .8, .995, .02}]
]

binList20[img_] := Module[{gimg = gradient3[img]},
  Table[DeleteSmallComponents[Erosion[Binarize[gimg, FindThreshold[img, Method -> {"BlackFraction", t}]], 1]], {t, .8, .995, .01}]
]

binList40[img_] := Module[{gimg = gradient3[img]},
  Table[DeleteSmallComponents[Erosion[Binarize[gimg, FindThreshold[img, Method -> {"BlackFraction", t}]], 1]], {t, .8, .9975, .005}]
]

morphologicalEulerNumbers10[img_] := MorphologicalEulerNumber /@ binList10[img]
morphologicalEulerNumbers20[img_] := MorphologicalEulerNumber /@ binList20[img]
morphologicalEulerNumbers40[img_] := MorphologicalEulerNumber /@ binList40[img]

gradient2Persistence[img_] := persistence[64, gradient2, img]
gradient3Persistence[img_] := persistence[64, gradient3, img]
gradient4Persistence[img_] := persistence[64, gradient4, img]

allFeatures[img_] := ImageAdjust /@ KarhunenLoeveDecomposition[ImageAdjust /@ Flatten[ColorSeparate /@ Flatten[{
  ColorConvert[img, "LUV"],
  ColorConvert[img, "RGB"],
  Table[GradientFilter[img, i], {i, 1, 5, 2}],
  Table[CornerFilter[img, i], {i, 2, 4, 1}],
  Table[RidgeFilter[img, i], {i, 2, 4, 1}],
  ImagePeriodogram[img],
  Table[LaplacianGaussianFilter[img, i], {i, 2, 4, 1}]
}]]][[1]]

kldCombined[img_] := ColorCombine[allFeatures[img][[;; 3]]]

kldGradients[img_] := ColorCombine[(ImageAdjust /@ KarhunenLoeveDecomposition[ColorConvert[#, "Grayscale"]& /@ Join[ColorSeparate[ColorConvert[img, "LUV"]], Table[ImageAdjust[GradientFilter[img, i]], {i, 1, 8, 3}]]][[1]])[[;; 3]]]

functions = SortBy[{Identity, luma, grayscale, gradient2, gradient4, gabor2x, gabor2y, gabor4x, gabor4y, corner2, corner4, ridge2, ridge4, laplacianGaussianFilter1, laplacianGaussianFilter2, laplacianGaussianFilter3, laplacian1, laplacian2, gradient2Persistence, gradient3Persistence, gradient4Persistence}, ToString];

methods = Sort[{
  "LogisticRegression",
  "NaiveBayes",
  {"NearestNeighbors", "NeighborsNumber" -> 25},
  {"NearestNeighbors", "NeighborsNumber" -> 50},
  {"NearestNeighbors", "NeighborsNumber" -> 100},
  {"RandomForest", "TreeNumber" -> 50, "LeafSize" -> 5},
  {"RandomForest", "TreeNumber" -> 100, "LeafSize" -> 10},
  {"RandomForest", "TreeNumber" -> 200, "LeafSize" -> 20}
}];