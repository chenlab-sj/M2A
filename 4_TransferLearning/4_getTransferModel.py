#!/usr/bin/python

'''
Requires:
python/3.6.5+
keras
tensorflow
h5py
scipy
os
pprint
pandas
numpy

Purpose:
Perform transfer learning on a single sample, updating a current model

Input Required:
	1) H5 feature file.
	2) Model weights

Output:
	1) Updated model weights
'''
import numpy as np
np.random.seed(9)
import os
import h5py
import pprint
import argparse
import pandas as pd
from scipy.stats import pearsonr
from sklearn.utils import shuffle
from keras.callbacks import EarlyStopping
from keras.models import Sequential
from keras.models import load_model
from keras.callbacks import Callback as ClassCallBack

#------------------------------------------------------------------------------
def parseArguments(): 
	# Create argument parser
	parser = argparse.ArgumentParser()   
	# Positional mandatory arguments
	parser.add_argument("FeatureFilePath", 
		help="Full path to feature .h5 file.", type=str)
	parser.add_argument("ModelFilePath", 
		help="Full path to model .h5 file.", type=str)
	parser.add_argument("--outFileName", help="File name prefix to use for output", type=str)
	parser.add_argument("--outDirectory", help="Directory to write output to", default="output", type=str)
	# Parse arguments
	args = parser.parse_args()
	return args

#------------------------------------------------------------------------------
# 2.1.1) build LossHistory object
class getLossHistory(ClassCallBack):
	def on_train_begin(self, logs={}):
		self.losses = []

	def on_epoch_end(self, epoch, logs={}):
		self.losses.append(logs.get('val_loss'))
	
	def early_stop_epoch(self):
		return [i for i,l in enumerate(self.losses) if l == min(self.losses)][0]

#------------------------------------------------------------------------------
# 2.1 build callback classes
def getCallBacks():
	# early stopping returns the best model performance after
	# 10 steps, providing no additional improvement (.0001)
	ES = EarlyStopping(monitor='val_loss',
						min_delta = .0001,
						verbose=0,
						mode="min",#3 "auto" would detect loss or acc
						patience=10,
						baseline=None,
						restore_best_weights=True)
	# 2.1.1) build LossHistory object
	LossHistory = getLossHistory()
	return ES, LossHistory
	
#------------------------------------------------------------------------------
# 1.1) parse H5 file for input to model
def parseInput(FeatureFilePath):
	SampleDataHDF = h5py.File(FeatureFilePath, mode="r")
	# train data
	TrainData = np.array(SampleDataHDF["FeatureInput"])
	# respvar
	TrainRespVar = np.array(SampleDataHDF["log2_ChipDivInput"])
	# shuffle data
	TrainData, TrainRespVar = shuffle(TrainData, 
									TrainRespVar,
									random_state=0)
	return TrainData, TrainRespVar

#------------------------------------------------------------------------------
# 1) Initalize transfer learning object
class learning_pipeline(object):
	#-------------------------------------------------------------------------#
	def __init__(self, FeatureFilePath, ModelFilePath, OutputFilePath):
		self.OutputFilePath = OutputFilePath
		self.ModelFilePath = ModelFilePath
		# 1.1) parse H5 file for input to model
		self.TrainData, self.TrainRespVar = parseInput(FeatureFilePath)
		# Training parameters
		self.Valid_Split = .2
		self.Epochs = 80
		self.Batch = 128
		self.Verb = 0
	#-------------------------------------------------------------------------#
	# 2) Training and return model instances
	def trainCNN(self):
		# load base model
		BaseModel = load_model(self.ModelFilePath)
		# freeze cnn layers from basemodel
		for Layer in BaseModel.layers[:-3]:
			Layer.name = "BaseModel_" + Layer.name 
			Layer.trainable = False
		# add basemodel layers and compile
		TransferModel = Sequential()
		TransferModel.add(BaseModel)
		TransferModel.compile(loss='mse', optimizer="adadelta",metrics=['mse'])
		# 2.1) build callback classes
		CallBacks, LossHistory = getCallBacks()
		# fitting model to unfrozen fc-nn layers
		FitModel = TransferModel.fit(self.TrainData, 
								self.TrainRespVar, 
								validation_split=self.Valid_Split, 
								epochs=self.Epochs, 
								batch_size=self.Batch,
								callbacks=[CallBacks, LossHistory], 
								verbose=self.Verb)			
		return FitModel, TransferModel, LossHistory

	#-------------------------------------------------------------------------#
	# 3) get training performance:
	def getTrainPerformance(self, TransferModel, FitModel):
		#---------------------------------------------------------------------#
		#3.1) get R2, RMSE
		def getStatistics(Pred, Actual, OutputFilePath):
			# get R^2, RMSE
			R, P = pearsonr(Pred, Actual)
			RMSE = np.sqrt(((Pred - Actual) ** 2).mean())
			# format extra dec places
			Stats = [format(v,".3f") for v in [R**2, P, RMSE]]
			# return formated stats dict text
			StatsDict = dict(R2=Stats[0],pval=P, RMSE=Stats[1])
			StatsText = pprint.pformat(StatsDict,width=20)
			StatsText = StatsText.replace("{","").replace("}","")
			print("Training Statistics:","\n", StatsText)
			# print to file
			FileName = OutputFilePath+"_TrainStatistics.txt"
			open(FileName, 'w').write(StatsText)
			return StatsDict

		#---------------------------------------------------------------------#
		# 3.2) plot the training/validation curve
		def plotRegrLossCurve(OutputFilePath, FitModel):
			import matplotlib
			matplotlib.use('Agg')
			import matplotlib.pyplot as plt
			#-----------------------------------------------------------------#
			PlotNamePath = OutputFilePath+"_Model_Loss_Curve.png"
			# plot loss
			Fig, Ax = plt.subplots(ncols=1,nrows=1, dpi=80, 
									facecolor='w', edgecolor='k')
			plt.plot(FitModel.history['loss'])
			plt.plot(FitModel.history['val_loss'])
			plt.title('Model Loss')
			plt.ylabel('Loss, MSE')
			plt.xlabel('Epoch')
			plt.legend(['train', 'validation'], loc='upper left')
			Fig.savefig(PlotNamePath)

		#---------------------------------------------------------------------#
		# return predictions on original training data:
		Predictions = TransferModel.predict(self.TrainData, 
										batch_size=250, 
										verbose=0)
		Predictions = Predictions.flatten()
		# print stats to stdout and to file
		getStatistics(Predictions, self.TrainRespVar, self.OutputFilePath)
		# save training curve to file
		plotRegrLossCurve(self.OutputFilePath, FitModel)

	#-------------------------------------------------------------------------#
	# 4) save transfer model
	def saveModel(self, FitModel, LossHistory):
		NumEpochs = str(LossHistory.early_stop_epoch())
		ModelPath = self.OutputFilePath+"_best_model_ep"+NumEpochs+".h5"
		FitModel.save(ModelPath)

#------------------------------------------------------------------------------
def main(LineArgs):
	# variable input
	FeatureFilePath = LineArgs.FeatureFilePath
	ModelFilePath = LineArgs.ModelFilePath
	# determine output directory
	Sample = ".".join(FeatureFilePath.split("/")[-1].split(".")[:-1])
	OutputFilePath = LineArgs.outDirectory
	if not os.path.exists(LineArgs.outDirectory):
		os.mkdir(LineArgs.outDirectory)
	# determine fullpath to output file path + prefix
	if not LineArgs.outFileName:
		OutputFilePath += "/" + "NewModel_"+Sample
	else:
		OutputFilePath += "/" + LineArgs.outFileName
	# 1) Initalize transfer learning object
	# 	1.1) parse H5 file for input to model
	TrainPipe = learning_pipeline(FeatureFilePath, ModelFilePath, OutputFilePath)
	# 2) Training and return model instances
	# 	2.1) build callback classes
	FitModel, TransferModel, LossHistory = TrainPipe.trainCNN()
	# 3) get training performance:
	# 	3.1) get R2, RMSE
	# 	3.2) plot the training/validation curve
	TrainPipe.getTrainPerformance(TransferModel, FitModel)
	# 4) save transfer model
	TrainPipe.saveModel(TransferModel, LossHistory)

# Capture command line args, with or without defaults
if __name__ == '__main__':
	# Parse the arguments
	LineArgs = parseArguments()
	main(LineArgs)