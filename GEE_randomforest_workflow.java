var model_df = ee.FeatureCollection("projects/phd-chapters/assets/sociable_lapwing/model_df_gee.csv"),
    BL_HR = ee.FeatureCollection("projects/phd-chapters/assets/BirdLife_HR"),
    roughness = ee.ImageCollection("projects/sat-io/open-datasets/Geomorpho90m/roughness"),
    regional_changes = ee.FeatureCollection("projects/grants-463915/assets/regional_changes");
	
//==========================================================================================================
// Setup
//==========================================================================================================
/*
- "model_df" is the file shared by Leon from HU where the values for NDVI and TC indices were 
   extracted for each presence and absence point. "obs_dat" column for that data is a date column which we dont 
   really need, so ignore that. The status column has the presence points along with 10 absence points
   that you generated. 
*/
//----------------------------------------------------------------------------------------------------------
// Function to add a human-readable date column based on the 'obs_dat' property.
//----------------------------------------------------------------------------------------------------------
var add_timestamp = function(feature) {
  // Set the formatted date as a new property.
  return feature.set('system:time_start', feature.get('obs_dat'));
};
model_df = model_df.map(add_timestamp); 
//==========================================================================================================
// define study region (roi)
//==========================================================================================================
var roi = BL_HR; // This is the BirdLife ranges for the species
Map.centerObject(roi, 6); 
//==========================================================================================================
// PREDICTION VARIABLES SETUP - LANDSAT BASED TCA + NDVI + ROUGHNESS
//==========================================================================================================
//==========================================================================================================
// Landsat archive data processing to generate predictor surface
//==========================================================================================================
// function to combine landsat imagery from landsat 4-8 into one collection
// and perform cloud masking based on pixel-qa masks 
function LandsatCollection(roi, start, end, masks){
  //-------------
  // Define masks
  masks = (typeof masks !== 'undefined') ?  masks : ee.List(['clouds', 'fill', 'water']);

  //----------------------------------------------------------
  // Collect, combine and process all available Landsat images
  var bands = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5','SR_B7','QA_PIXEL'];
  var band_names = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2','QA_PIXEL'];

  //----------------------------
  // For Landsat 8 and Landsat 9
  var l8_9_bands = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7','QA_PIXEL'];
  var l8_9_band_names = ['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2','QA_PIXEL'];
  
  //--------
  // QA bits 
  var mask_clouds = ee.List(masks).contains("clouds");
  var mask_water = ee.List(masks).contains("water");
  //var mask_snow = ee.List(masks).contains("snow");
  var mask_fill = ee.List(masks).contains("fill");
  
  var cloudbit = ee.Number(ee.Algorithms.If(mask_clouds, 40, 0));
  var waterbit = ee.Number(ee.Algorithms.If(mask_water, 4, 0));
  //var snowbit = ee.Number(ee.Algorithms.If(mask_snow, 16, 0));
  var fillbit = ee.Number(ee.Algorithms.If(mask_fill, 1, 0));
  
  //---------------------
  // Combine all together
  var bits = cloudbit.add(waterbit).add(fillbit)//.add(snowbit);
  
  //-----------------
  // Helper functions
  // function to apply masks based on pixel qa band
  var apply_masks = function(img){
    var qa = img.select('QA_PIXEL');
    var mask = qa.bitwiseAnd(bits).eq(0);
    return img.updateMask(mask);
  };

  //-----------------------------------
  // Apply functions on Landsat archive (4-9)
  var l4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2")
                                .select(bands, band_names)
                                .filterBounds(roi)
                                .filterDate(start, end);
                                
  var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2")
                                .select(bands, band_names)
                                .filterBounds(roi)
                                .filterDate(start, end);
                                
  var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2")
                                .select(bands, band_names)
                                .filterBounds(roi)
                                .filterDate(start, end);
                                
  var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                                .select(l8_9_bands, l8_9_band_names)
                                .filterBounds(roi)
                                .filterDate(start, end);
  var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2")
                                .select(l8_9_bands, l8_9_band_names)
                                .filterBounds(roi)
                                .filterDate(start, end);
                         
  //---------------------------------------------------------------------------------------------------------------------------------
  // Combine landsat collection
  var landsat = ee.ImageCollection(l4.merge(l5).merge(l7).merge(l8).merge(l9)).map(apply_masks);
  return landsat;
} // End of funtion LandsatCollection

//========================================================================================================================
// Tasseled cap Function + NDVI + Roughness
//========================================================================================================================
// function to calculate tasseled cap indices from collection of 
// landsat surface reflectance imagery
function IndicesCollection(LandsatCollection){
  var tc_ndvi = function(image) {
    var img = image.select(['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']);
    
    // coefficients for Landsat surface reflectance (Crist 1985)
    // Crist, E.P. (1985). A TM Tasseled Cap equivalent transformation for reflectance factor data.
    var brightness_c= ee.Image([0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303]);
    var greenness_c= ee.Image([-0.1603, 0.2819, -0.4934, 0.7940, -0.0002, -0.1446]);
    var wetness_c= ee.Image([0.0315,  0.2021,  0.3102,  0.1594, -0.6806, -0.6109]);
    
    var brightness = img.multiply(brightness_c);
    var greenness = img.multiply(greenness_c);
    var wetness = img.multiply(wetness_c);
    
    brightness = brightness.reduce(ee.call("Reducer.sum"));
    greenness = greenness.reduce(ee.call("Reducer.sum"));
    wetness = wetness.reduce(ee.call("Reducer.sum"));
    //---------------------------------------------------------------------------------------------------------------------------------
    // Add NDVI
    var ndvi = image.normalizedDifference(['NIR', 'R'])
                    .multiply(10000)
                    .rename('NDVI');
    //---------------------------------------------------------------------------------------------------------------------------------  
    // Add Roughness
    var roughness_mosaic = roughness.mosaic().rename('Roughness');
    //---------------------------------------------------------------------------------------------------------------------------------
    // Put all together
    var indices = ee.Image(brightness)
                      .addBands(greenness)
                      .addBands(wetness)
                      .addBands(ndvi)
                      .addBands(roughness_mosaic)
                      .rename(['TCB','TCG','TCW','NDVI','Roughness']).clip(BL_HR);
    
    return indices;
  };
  var out = LandsatCollection.map(tc_ndvi); // Plus roughness, although not in the name
  return out;
}

//==========================================================================================================
// Calculate spatial-temporal metrics
/* 
  Training data were summarized over three months (See R scripts) using median of TC values
  The predictor surfaces will follow the same time-intervals. We are only working with median values
*/
//==========================================================================================================
function calculate_metrics(roi, start, end) {

  // collect all landsat data in target period
  var landsat_data = LandsatCollection(roi, start, end, ['clouds','fill','water']); 

  // compute indices from collected data
  var tc = IndicesCollection(landsat_data);
  var out = tc.reduce(ee.Reducer.median())
              .rename(['TCB','TCG','TCW','NDVI','Roughness']);
  return out;
}

//---------------------------------------------------------------------------------------------------------------------------------
// indicate which years should covered in the monitoring
//---------------------------------------------------------------------------------------------------------------------------------
var years = ee.List.sequence(1995,2024); // NOT TO EXTRACT ANYTHING IN RASTER FORMAT AT THIS SETTING

// Create a list of all months for all years
var months = years.map(function(year) {
  return ee.List.sequence(1, 12).map(function(month) {
    return ee.Date.fromYMD(ee.Number(year), ee.Number(month), 1);
  });
}).flatten();

//---------------------------------------------------------------------------------------------------------
/* Calculate spectral-temporal metrics for EACH MONTH
   This is an important step. Here I choose month and then buffer it such that its a running window 
   of 3-month period. This closes some data gaps. More importantly, it matches the 
   time-calibration of the training data which is important. This way the prediction surfaces have the index
   values summarized in the same way the training data values are.
   This a bit weird here so read the code carefully. For training data we actually use +/- 45 days to make three 
   months. Here, we do +/- 30 because 30 days are counted before the FIRST day of the month and LAST day of the month
   which excludes the running month instead, and spanning roughly over 3-month window. 
*/
//---------------------------------------------------------------------------------------------------------
var monthly_list = months.map(function(start) {
  start = ee.Date(start); // Ensure start is an ee.Date object
  var end = start.advance(1, 'month'); // Define the end date
  
  // Call the calculate_metrics function
  return calculate_metrics(roi, start.advance(-30, 'day'), end.advance(30, 'day')) // Add buffer if needed
    .set('system:time_start', start.millis()) // Set a time property
    .set('month', start.format('YYYY-MM')) // Set a month identifier
    .set('year', start.format('YYYY')) // Set a year identifier 
    //.set('date', start.format('YYYY-MM-dd'))
});

// Convert the list of results into an ImageCollection
var monthly_metrics = ee.ImageCollection(monthly_list);

//---------------------------------------------------------------------------------------------------------
// Function to remove NA pixels uniformly from all the bands
//---------------------------------------------------------------------------------------------------------
var RemoveNA = function(image) {
  // Retrieve the mask for all bands and compute the minimum mask
  var combinedMask = image.mask().reduce(ee.Reducer.min());
  // Apply the combined mask to the image
  return image.updateMask(combinedMask);
};
monthly_metrics = monthly_metrics.map(RemoveNA);

// Print or display the collection
//print('Monthly Metrics Collection:', monthly_metrics); // THIS IS YOUR GLOBAL PREDICTION SURFACE
//Map.addLayer(monthly_metrics.toBands().clip(geometry), {}, "test");

//==========================================================================================================
// TRAINING DATA SETUP
//==========================================================================================================
//==========================================================================================================
// Adding Roughness to TRAINING DATA, because this was not extracted together with 
// TCA and NDVI
//==========================================================================================================
var roughness_mosaic = roughness.mosaic().rename('Roughness');

// Define a function to add roughness
var addValues = function(feature) {
  
  // Extract roughness value
  var roughnessValue = roughness_mosaic.reduceRegion({
    reducer: ee.Reducer.first(),
    geometry: feature.geometry(),
    scale: 30, // Match the dataset resolution
    maxPixels: 1e8
  }).get('Roughness'); // Adjust the band name if different
  
  // Add values as properties
  return feature.set({
    Roughness: roughnessValue,
  });
};

//----------------------------------------------------------------------------------------------------------
// Add Roughness to dataframe
//----------------------------------------------------------------------------------------------------------
// Map the function over the FeatureCollection to add the roughness column
var model_df = model_df.map(addValues); //print(model_df.limit(10), 'model_df_sample');

// Remove NAs 
model_df = model_df.filter(ee.Filter.and(ee.Filter.notNull(['Roughness'])));

//==========================================================================================================
// Model training with k-fold (10) cross-validation 
// (Based on model_df imported into GEE. NO IMAGERY INVOLVED HERE. Only model training)
//==========================================================================================================
//----------------------------------------------------------------------------------------------------------
// Define predictors and response variables
//----------------------------------------------------------------------------------------------------------
var predictors = ['NDVI','TCB','TCG','TCW','Roughness']; // Feature names
var response = 'status'; // Response variable

//----------------------------------------------------------------------------------------------------------
// Function to split a feature collection into k-folds
//----------------------------------------------------------------------------------------------------------
var kFoldSplitStratified = function(features, folds) {
    var step = ee.Number(1).divide(folds);
    var thresholds = ee.List.sequence(
        0, ee.Number(1).subtract(step), step);

    // Separate the features into presence (1) and absence (0) classes
    var presence = features.filter(ee.Filter.eq('status', 1));  // Assume 'class' is the property name for presence/absence
    var absence = features.filter(ee.Filter.eq('status', 0));

    // Shuffle both the presence and absence features
    presence = presence.randomColumn({seed: 0});
    absence = absence.randomColumn({seed: 0});

        // Perform the stratified split by applying the same thresholds to both classes
    var splits = thresholds.map(function (threshold) {
        var trainingPresence = presence.filter(
            ee.Filter.or(
                ee.Filter.lt('random', threshold),
                ee.Filter.gte('random', ee.Number(threshold).add(step))
            )
        );
        var validationPresence = presence.filter(
            ee.Filter.and(
                ee.Filter.gte('random', threshold),
                ee.Filter.lt('random', ee.Number(threshold).add(step))
            )
        );
        var trainingAbsence = absence.filter(
            ee.Filter.or(
                ee.Filter.lt('random', threshold),
                ee.Filter.gte('random', ee.Number(threshold).add(step))
            )
        );
        var validationAbsence = absence.filter(
            ee.Filter.and(
                ee.Filter.gte('random', threshold),
                ee.Filter.lt('random', ee.Number(threshold).add(step))
            )
        );

        // Combine the presence and absence splits
        var trainingSplit = trainingPresence.merge(trainingAbsence);
        var validationSplit = validationPresence.merge(validationAbsence);

        return ee.Feature(null, {
            'training': trainingSplit,
            'validation': validationSplit
        });
    });

    return ee.FeatureCollection(splits);
};

// Use k=10
var folds = kFoldSplitStratified(model_df, 5);
print(folds);

// Check the split
//print('Fold 0 Training Samples', folds.first().get('training'));
//print('Fold 0 Validation Samples', folds.first().get('validation'));
//===============================================================================================================
// Accuracy assessment
//===============================================================================================================
// Assess the accuracy for each pair of training and validation
// THe "folds" are taken from the k-fold cross validation step above, so they are connected
/*
//---------------------------------------------------------------------------------------------------------------
// This is default written in the original workflow
//---------------------------------------------------------------------------------------------------------------
var accuracies = ee.FeatureCollection(folds.map(function(fold) {
  var trainingGcp = ee.FeatureCollection(fold.get('training'));
  var validationGcp = ee.FeatureCollection(fold.get('validation'));

  // Train a classifier.
  var classifier = ee.Classifier.smileRandomForest({
    numberOfTrees: 150,
    variablesPerSplit: null, 
    minLeafPopulation: 10,
    maxNodes: null
  }).setOutputMode('PROBABILITY') // Output probabilities
    .train({
      features: trainingGcp,  
      classProperty: 'status',
      inputProperties: predictors
    });
  
  // Validate the classifier.
  var validation = validationGcp.classify(classifier); 
  
  // Threshold probabilities to classify them as 0 or 1
  var threshold = 0.5; // Default threshold
  var classifiedValidation = validation.map(function(feature) {
    var probability = feature.get('classification'); // Probability of being 1
    var predictedClass = ee.Number(probability).gt(threshold) // Compare with threshold
      .int(); // Convert to integer (0 or 1)
    return feature.set('predicted_class', predictedClass);
  });

  // Compute the accuracy using the predicted_class field.
  var accuracy = classifiedValidation
    .errorMatrix('status', 'predicted_class')
    .accuracy();

  return ee.Feature(null, {'accuracy': accuracy});
}));

// Print and export results
// print('K-fold Validation Results', accuracies.aggregate_array('accuracy'));
print('K-fold Validation Results', accuracies);

var exportTable = accuracies.map(function(feature) {
  return ee.Feature(null, {
    'accuracy': feature.get('accuracy')
  });
});

Export.table.toDrive({
  collection: exportTable,
  description: 'KFoldValidationAccuracies',
  fileFormat: 'CSV'
});
*/
//---------------------------------------------------------------------------------------------------------------
// Tweaked above part using GEE to export predicted values from all k-folds. ROC-AOC is then calculated in R using 
// pROC and ggplot2 packages
//---------------------------------------------------------------------------------------------------------------
var allFoldResults = ee.FeatureCollection(folds.map(function(fold) {
  var trainingGcp = ee.FeatureCollection(fold.get('training'));
  var validationGcp = ee.FeatureCollection(fold.get('validation'));

  var classifier = ee.Classifier.smileRandomForest({
    numberOfTrees: 200,
    variablesPerSplit: null,
    minLeafPopulation: 3
  }).setOutputMode('PROBABILITY')
    .train({
      features: trainingGcp,
      classProperty: 'status',
      inputProperties: predictors
    });

  // Classify the validation set to get probabilities
  var validated = validationGcp.classify(classifier);

  // Extract probability and true label
  return validated.map(function(f) {
    return ee.Feature(null, {
      'status': f.get('status'),
      'probability': f.get('classification')
    });
  });
}).flatten()); // Flatten to merge all folds together

// print('All Fold Results Size:', allFoldResults);

Export.table.toDrive({
  collection: allFoldResults,
  description: 'KFold_RF_Validation_ProbabilitiesLeaf3_Trees1000',
  fileFormat: 'CSV'
});

//---------------------------------------------------------------------------------------------------------------
// Adapted from Ramiro D.Crego SDM-GEE paper
//---------------------------------------------------------------------------------------------------------------
// Define functions to estimate sensitivity, specificity and precision.
var getROCMetrics = function(resultsFC) {
  var thresholds = ee.List.sequence(0, 1, 0.04); // 25 thresholds
  return ee.FeatureCollection(thresholds.map(function(thresh) {
    thresh = ee.Number(thresh);

    var pres = resultsFC.filter(ee.Filter.eq('status', 1));
    var abs = resultsFC.filter(ee.Filter.eq('status', 0));

    var TP = pres.filter(ee.Filter.gte('probability', thresh)).size();
    var FN = pres.filter(ee.Filter.lt('probability', thresh)).size();
    var FP = abs.filter(ee.Filter.gte('probability', thresh)).size();
    var TN = abs.filter(ee.Filter.lt('probability', thresh)).size();

    var TPR = ee.Number(TP).divide(ee.Number(TP).add(FN)); // Sensitivity
    var FPR = ee.Number(FP).divide(ee.Number(FP).add(TN)); // 1 - Specificity

    return ee.Feature(null, {
      'threshold': thresh,
      'TPR': TPR,
      'FPR': FPR
    });
  }));
};

var computeAUC = function(rocFC) {
  var FPR = ee.Array(rocFC.aggregate_array('FPR'));
  var TPR = ee.Array(rocFC.aggregate_array('TPR'));

  var dFPR = FPR.slice(0, 1).subtract(FPR.slice(0, 0, -1));
  var avgTPR = TPR.slice(0, 1).add(TPR.slice(0, 0, -1)).multiply(0.5);

  return dFPR.multiply(avgTPR).reduce('sum', [0]).abs().get([0]);
};

var roc = getROCMetrics(allFoldResults);
// print('ROC-AUC (training data only):', roc);
var auc = computeAUC(roc);
print('ROC-AUC (training data only):', auc);

// Export ROC data as CSV
Export.table.toDrive({
  collection: roc,
  description: 'ROC_Curve',
  fileFormat: 'CSV'
});


//===============================================================================================================
// Hyperparameter tuning
//===============================================================================================================
// Tune the numberOfTrees parameter.

var numTreesList = ee.List.sequence(50, 200, 10); //Try from 10 - 200, with 30 tree increment
var withRandom = model_df.randomColumn('random'); // Add random column
var trainingSet = withRandom.filter(ee.Filter.lt('random', 0.7)); // 70% for training
var validationSet = withRandom.filter(ee.Filter.gte('random', 0.7)); // 30% for validation

var accuracies = numTreesList.map(function(numTrees) {
  var classifier = ee.Classifier.smileRandomForest(numTrees)//.setOutputMode('PROBABILITY')
      .train({
        features: trainingSet,
        classProperty: 'status',
        inputProperties: predictors
      });

  // Here we are classifying a table instead of an image
  // Classifiers work on both images and tables
  return validationSet.classify(classifier)
      .errorMatrix('status', 'classification')
      .accuracy();
});


var chart = ui.Chart.array.values({
  array: ee.Array(accuracies),
  axis: 0,
  xLabels: numTreesList
  }).setOptions({
      title: 'Hyperparameter Tuning for the numberOfTrees Parameters',
      vAxis: {title: 'Validation Accuracy'},
      hAxis: {title: 'Number of Tress', gridlines: {count: 15}}
  });

print(chart);
//===============================================================================================================
// Model prediction (Here we need to start filtering monthly metrics, if we want to predict on seasonal surfaces)
//===============================================================================================================
/*
  - Predictions are still made on monthly level 
  - They are summarized AFTER the habitat suitability probability is acquired per month
*/
//===============================================================================================================
// Apply the classifier to each image in the collection
// Train a classifier.
/*
var classifier_df = ee.Classifier.smileRandomForest({
    numberOfTrees: 150,
    variablesPerSplit: null, 
    minLeafPopulation: 10,
    maxNodes: null
  }).setOutputMode('PROBABILITY') // Output probabilities
    .train({
      features: trainingSet,  
      classProperty: 'status',
      inputProperties: predictors
    });
    
// All seasons Monthly - 400+ maps are generated here. We map here BECAUSE MONTHLY METRICS IS AN IMAGE COLLECTION
var classifiedCollection_df = monthly_metrics.map(function(image) {
  return image.classify(classifier_df) 
              .set('month', image.get('month')) // Optional: retain month metadata
              .set('year', image.get('year'))
              .set('system:time_start', image.date().millis()); // Must have date info in this format 
});

//===============================================================================================================
// Exporting results for ROC-AUC calculation (Initial attempt at ROC-AUC. Not implemented atm)
//===============================================================================================================
// // Classify the validation set to get probabilities
//var validated_df = validationSet.classify(classifier_df, 'predicted_prob');

// Export as CSV
//Export.table.toDrive({
//  collection: validated_df,
//  description: 'RF_Probabilities_ValidationSet',
//  fileFormat: 'CSV'
//});
//=================================================================================================
// Feature Importance
//=================================================================================================
// Run .explain() to see what the classifer looks like
//print(classifier.explain());
// Calculate variable importance
/*
var importance = ee.Dictionary(classifier_df.explain().get('importance'))

// Calculate relative importance
var sum = importance.values().reduce(ee.Reducer.sum())

var relativeImportance = importance.map(function(key, val) {
   return (ee.Number(val).multiply(100)).divide(sum)
  })
//print(relativeImportance)

// Create a FeatureCollection so we can chart it
var importanceFc = ee.FeatureCollection([
  ee.Feature(null, relativeImportance)
])

var chart = ui.Chart.feature.byProperty({
  features: importanceFc
}).setOptions({
      title: 'Feature Importance',
      vAxis: {title: 'Importance'},
      hAxis: {title: 'Feature'}
  })
//print(chart)
/*
//----------------
// Breeding season
var classifiedCollection_Br = monthly_metrics.map(function(image) {
  return image.classify(classifier_Br)
              .set('month', image.get('month')) // Optional: retain month metadata
              .set('year', image.get('year'))
              .set('system:time_start', image.date().millis()); // Must have date info in this format 
});
//--------------
// Winter season
var classifiedCollection_Wn = monthly_metrics.map(function(image) {
  return image.classify(classifier_Wn)
              .set('month', image.get('month')) // Optional: retain month metadata
              .set('year', image.get('year'))
              .set('system:time_start', image.date().millis()); // Must have date info in this format 
});
//-------------------
// Migratory movement
var classifiedCollection_Mg = monthly_metrics.map(function(image) {
  return image.classify(classifier_Mg)
              .set('month', image.get('month')) // Optional: retain month metadata
              .set('year', image.get('year'))
              .set('system:time_start', image.date().millis()); // Must have date info in this format 
});
*/
// Print to inspect
// print('Classified Collection', classifiedCollection_df);
// print('Classified Collection:Breeding', classifiedCollection_Br);

//==========================================================================================================
// Summarizing monthly probabilities to interpret habitat suitability through time
   //-----------------------------------------------------------------------------------
/* Output 1 - Decedal (1995-2024) Net change in annual habitat suitability (using the entire BirdLife range map)
   Based on monthly output (classifiedCollection_df) decedal change in habitat suitability is compared using 
   6 yearly median Habitat suitability maps representing net change in suitability for that decade.
   //-----------------------------------------------------------------------------------
   Output 2 - Regional trends and changes in habitat suitability are calculated for the 15 test regions representing 
   Breeding, stop-over and wintering home range of the Sociable lapwing. In this case yearly medians are extracted for the
   regional 100km² areas. They are ploted using locally estimated scatter plot smoothing (loess) to visualize 
   increasing and decreasing patterns in habitat suitability
*/
//==========================================================================================================
//----------------------------------------------------------------------------------------------------------
// Output 1
//----------------------------------------------------------------------------------------------------------
var visParams = {
  min: 0,
  max: 0.4,
  palette: ['yellow','green'],
};
//----------------------------------------------------------------------------------------------------------
// 1995-2004
var caf_1995_median = classifiedCollection_df.filter(ee.Filter.date('1994-01-01','1996-12-31')).median().rename('HSI_median_1995'); 
var caf_2004_median = classifiedCollection_df.filter(ee.Filter.date('2004-01-01','2004-12-31')).median().rename('HSI_median_2004'); 

//-----------
// 2005-2014
var caf_2005_median = classifiedCollection_df.filter(ee.Filter.date('2005-01-01','2005-12-31')).median().rename('HSI_median_2005'); 
var caf_2014_median = classifiedCollection_df.filter(ee.Filter.date('2014-01-01','2014-12-31')).median().rename('HSI_median_2014'); 

//----------
// 2015-2024
var caf_2015_median = classifiedCollection_df.filter(ee.Filter.date('2015-01-01','2015-12-31')).median().rename('HSI_median_2015');
var caf_2024_median = classifiedCollection_df.filter(ee.Filter.date('2024-01-01','2024-12-31')).median().rename('HSI_median_2024');

//----------
// 1995-2004
Export.image.toDrive({
  image: caf_1995_median,
  description: 'HSI_1995',
  scale: 1000,
  region: BL_HR, // Define your region
  maxPixels: 1e13
});
Export.image.toDrive({
  image: caf_2004_median,
  description: 'HSI_2004',
  scale: 1000,
  region: BL_HR, // Define your region
  maxPixels: 1e13
});

//----------
// 2005-2014
Export.image.toDrive({
  image: caf_2005_median,
  description: 'HSI_2005',
  scale: 1000,
  region: BL_HR, // Define your region
  maxPixels: 1e13
});
Export.image.toDrive({
  image: caf_2014_median,
  description: 'HSI_2014',
  scale: 1000,
  region: BL_HR, // Define your region
  maxPixels: 1e13
});

//----------
// 2015-2024
Export.image.toDrive({
  image: caf_2015_median,
  description: 'HSI_2015',
  scale: 1000,
  region: BL_HR, // Define your region
  maxPixels: 1e13
});
Export.image.toDrive({
  image: caf_2024_median,
  description: 'HSI_2024',
  scale: 1000,
  region: BL_HR, // Define your region
  maxPixels: 1e13
});

//----------------------------------------------------------------------------------------------------------
// Output - 2
//----------------------------------------------------------------------------------------------------------
// Function to compute median per year
Map.addLayer(regional_changes, {}, "regional_changes")
var yearlyMedianCollection = ee.ImageCollection(
  years.map(function(year) {
    year = ee.Number(year);
    
    // Filter the monthly collection for the current year
    var yearlyImages = classifiedCollection_df.filter(ee.Filter.calendarRange(year, year, 'year'));
    
    // Take median across all images in the year
    var medianImage = yearlyImages.median()
      .set('year', year)
      .set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis());
    
    return medianImage;
  })
);

//----------------------------------------------------------------------------------------------------------
// Calculate regional median for each year betwen 1990 and 2024
//----------------------------------------------------------------------------------------------------------
var regional_changes_stats = yearlyMedianCollection.map(function(img) {
  var year = img.get('year');
  
  var stats = img.reduceRegions({
    collection: regional_changes,
    reducer: ee.Reducer.median(),
    scale: 1000
  });

  // Attach the year to each feature
  stats = stats.map(function(f) {
    return f.set('year', year);
  });

  return stats;
}).flatten();

//----------------------------------------------------------------------------------------------------------
// Export the file 
//----------------------------------------------------------------------------------------------------------
Export.table.toDrive({
  collection: regional_changes_stats,
  description: 'regional_changes_stats_V',
  fileFormat: 'CSV'
});

