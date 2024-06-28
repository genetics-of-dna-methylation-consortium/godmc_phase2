from collections import Counter

import numpy as np
import os
import tables
import gc
import pandas as pd
import Cryptodome.Random.random as crandom

class Encoder(object):

	def __init__(self, out):
		self.name='Encoder'
		self.N_matrix=None
		self.F=None
		self.F_inv=None
		self.encoder_chunks={}
		self.hdf5_iter=0
		self.npy_iter=Counter()
		self.pytable_filters = tables.Filters(complevel=9, complib='zlib')
		self.out=out
		self.metadata=None
		self.study_name=None
		self.phen_info_dic={}
		self.phen_info_dic['id']=None
		try:
			print ('Creating directories...')
			os.mkdir(os.path.join(self.out,'encode_genotype'))
			os.mkdir(os.path.join(self.out,'encode_phenotype'))
			os.mkdir(os.path.join(self.out,'encode_individuals'))
			os.mkdir(os.path.join(self.out,'encode_interaction'))

		except:
			print('Directories "encode_genotype","encode_phenotype","encode_individuals" are already exist in {}...'.format(self.out))


	def matrix(self,N, save=False):

		# Use cryptographic level randomness made available through the PyCryptodome Python package,
		# instead of the numpy.random library, which is completely deterministic.
		self.F = np.matrix(
			[crandom.randint(1, 9) for i in range(N * N)]).reshape(N, N)
		# Use 9 instead of 10 as the upper bound as this value is inclusive for PyCryptodome,
		# and exclusive for numpy.random.randint().

		# Old code for generating self.F. Now commented out since we use pycryptodome:
		#self.F=np.random.randint(1,10,N*N).reshape(N,N)
		self.F_inv=np.linalg.inv(self.F)
		self.F=np.array(self.F,dtype=np.float)
		if save:
			print ('Saving decoders Matrices to...{}'.format(self.out))
			np.save(os.path.join(self.out,'F.npy'),self.F)
			np.save(os.path.join(self.out,'F_inv.npy'),self.F_inv)



	def encode(self,data, data_type=None):

		if isinstance(self.F,type(None)):
			raise ValueError('Encode matrix is not define')

		x,y=data.shape
		#TODO (mid) add to encode_chunk info about data size and/or names indexes
		#TODO (mid) check the protocol if exist... all study on single stage

		a,b=self.F.shape

		if isinstance(data_type,type(None)):
			raise ValueError('data_type is None')
		elif data_type=="genotype":
			# If the data is genotype data, encode using F
			print ('Encoding genotype data')
			return np.dot(data,self.F)
		elif data_type=='phenotype':
			# If the data is phenotype data, encode using the inverse of F.
			print ('Encoding phenotype data')
			return np.dot(self.F_inv,data)


	def save_npy(self,data, save_path=None, info=None, index=None):
		if isinstance(save_path,type(None)):
			raise ValueError('Path cannot be none {}'.format(save_path))
		save_path = os.path.join(self.out, save_path)
		if os.path.isdir(save_path) is False:
			os.mkdir(save_path)
		np.save(os.path.join(save_path,str(self.npy_iter[save_path])+'_'+self.study_name+ '.npy'),data )
		if self.phen_info_dic['id'] is None:
			self.phen_info_dic['id']=np.array(info._data.id)[index]
		self.phen_info_dic[str(self.npy_iter[save_path])+'_'+self.study_name+ '.npy']=info._data.names[info._data.start:info._data.finish]
		self.npy_iter[save_path] += 1

	def save_csv(self,data,save_path=None, info=None, index=None):
		if isinstance(save_path,type(None)):
			raise ValueError('Path cannot be none {}'.format(save_path))
		save_path = os.path.join(self.out, save_path)
		if os.path.isdir(save_path) is False:
			os.mkdir(save_path)
		df=pd.DataFrame(data)
		df.columns=info._data.names[info._data.start:info._data.finish]
		df.insert(0,'id',np.array(info._data.id)[index])
		df.to_csv(os.path.join(save_path,str(self.npy_iter[save_path])+'_'+self.study_name+ '.csv'), sep='\t', index=False)
		self.npy_iter[save_path] += 1



	def save_hdf5(self,data, save_path=None,info=None, index=None):
		#only for genetics data

		if isinstance(save_path,type(None)):
			raise ValueError('Path cannot be none {}'.format(save_path))
		save_path = os.path.join(self.out, save_path)
		if os.path.isdir(save_path) is False:
			os.mkdir(save_path)

		print ('Saving data to ... {}'.format(os.path.join(save_path,str(self.hdf5_iter)+'_'+self.study_name+'_encoded.h5')))

		if not os.path.isfile(os.path.join(self.out,'encode_individuals',self.study_name + '.h5'  )):
			chunk = pd.DataFrame.from_dict({"individual": info._data.id[index] })
			chunk.to_hdf(os.path.join(self.out,'encode_individuals',self.study_name + '.h5' ), key='individuals', format='table', min_itemsize=25, complib='zlib', complevel=9)

		if isinstance(save_path,type(None)) or  not os.path.isdir(save_path):
			raise ValueError('There is no such path or directory {}'.format(save_path))

		self.h5_gen_file = tables.open_file(
			os.path.join(save_path,str(self.hdf5_iter)+'_'+self.study_name +'.h5'), 'w', title='encode_genotype')
		self.hdf5_iter+=1

		atom = tables.Float64Atom()
		self.genotype = self.h5_gen_file.create_carray(self.h5_gen_file.root, 'genotype', atom,
													  (data.shape),
													  title='Genotype', filters=self.pytable_filters)
		self.genotype[:] = data
		self.h5_gen_file.close()
		gc.collect()



