
from cnn_model import build_model, train_model
from .input import read_input_file
from .molecule_tensor import get_molecule_tensor
import os
import rmgpy
import numpy as np
from .data import get_HC_polycyclics_data_from_db, prepare_folded_data, prepare_data_one_fold
import logging

class Predictor(object):

	def __init__(self, input_file=None):

		self.model = None
		if input_file:
			self.input_file = input_file
		else:
			self.input_file = os.path.join(os.path.dirname(rmgpy.__file__),
										'cnn_framework',
										'data',  
										'predictor_input.py'
										)

	def build_model(self):
		"""
		This method is intended to provide a way to build default model 
		"""

		self.model = build_model()

	def load_input(self, path=None):
		"""
		This method is intended to provide a way to build model from an input file
		"""
		
		if path is None: 
			path = self.input_file
			print path
		read_input_file(path, self)

	def kfcv_train(self, folds):

		# prepare data for training
		(X, y) = get_HC_polycyclics_data_from_db('sdata134k', 'sdata134k_table')
		(folded_Xs, folded_ys) = prepare_folded_data(X, y, folds)

		losses = []
		val_losses = []
		for fold in range(folds):
			data = prepare_data_one_fold(folded_Xs, folded_ys, current_fold=fold)

			# execute train_model
			(model, loss, val_loss) = train_model(self.model, data, 
												nb_epoch=150, lr_func='0.01', 
												patience=10)

			losses.append(loss)
			val_losses.append(val_loss)
			# save model
			# self.save_model()

		# mean loss and val_loss used for selecting parameters, 
		# e.g., lr, epoch, attributes, etc
		mean_loss = mean(losses)
		mean_val_loss = mean(val_losses)

		self.report_training_results(mean_loss, mean_val_loss)

	def load_parameters(self, param_path=None):

		if not param_path:
			param_path = os.path.join(os.path.dirname(rmgpy.__file__),
									'cnn_framework',
									'data', 
									'weights', 
									'polycyclic_enthalpy_weights.h5'
									)

		self.model.load_weights(param_path)

	def save_parameters(self):

		pass

	def predict(self, molecule):

		molecule_tensor = get_molecule_tensor(molecule)

		molecule_tensor_array = np.array([molecule_tensor])
		return self.model.predict(molecule_tensor_array)[0][0]
    
