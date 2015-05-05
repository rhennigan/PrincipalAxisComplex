(* Mathematica source file  *)
(* Created by IntelliJ IDEA *)
(* :Author: rhennigan *)
(* :Date: 5/4/2015 *)

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