"""

Offline Standalone version of PDBMD2CD

Attributes:
    dsspfile_refs (list): Description
    pcdfiles_refs (list): Description
    ref_basis (list): Description
    ref_cd (list): Description
    ref_lincomb_ss (list): Description
    ref_lstsq_ss (list): Description

"""
import re
import string
import os
import time
import numpy as np
import sys
import math
import random
from math import sqrt
from scipy.spatial.distance import cdist, euclidean, pdist, squareform
from scipy.optimize import lsq_linear, minimize
from itertools import groupby, repeat

import argparse
#import requests
import time
import datetime

import tarfile
import zipfile
from shutil import make_archive, move, copy

import networkx as nx
import json
from collections import defaultdict

from Bio.PDB import PDBList, PDBIO, PDBParser, MMCIFParser

from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import subprocess

import config

twistcutoff=[-10,59]
an=[35, 110]
no_cand=12




def extract_file(path, to_directory='.'):
	"""Summary
	
	Args:
	    path (TYPE): Description
	    to_directory (str, optional): Description
	
	Raises:
	    ValueError: Description
	"""
	if path.endswith('.zip'):
		opener, mode = zipfile.ZipFile, 'r'
	elif path.endswith('.tar.gz') or path.endswith('.tgz'):
		opener, mode = tarfile.open, 'r:gz'
	elif path.endswith('.tar.bz2') or path.endswith('.tbz'):
		opener, mode = tarfile.open, 'r:bz2'
	else: 
		raise ValueError("Could not extract `%s` as no appropriate extractor is found" % path)
	
	cwd = os.getcwd()
	os.chdir(to_directory)
	
	try:
		file = opener(path, mode)
		try: file.extractall()
		finally: file.close()
	finally:
		os.chdir(cwd)


def pdbfile_to_dssp_local(directory, filename):
	"""Summary
	
	Args:
	    directory (TYPE): Description
	    filename (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	subprocess.call([config.DSSPPATH, "-i", os.path.join(directory,str(filename)), "-o", os.path.join(directory,str(filename)+".dssp")])

	try:
		with open(os.path.join(directory, str(filename)+".dssp"), "r") as dssp_in:
			return(dssp_in.readlines(), True)
	except FileNotFoundError:
		failed.append(filename)
		return("", False)



def angleTwoLines(a1, b1, c1, a2, b2, c2):
	"""Summary
	
	Args:
	    a1 (TYPE): Description
	    b1 (TYPE): Description
	    c1 (TYPE): Description
	    a2 (TYPE): Description
	    b2 (TYPE): Description
	    c2 (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	numerator=a1*a2+b1*b2+c1*c2
	denom=sqrt(a1**2+b1**2+c1**2)*sqrt(a2**2+b2**2+c2**2)
	costheta=numerator/denom
	angle=np.arccos(costheta)
	return(angle)

class DSSPData:

	"""Summary
	
	Attributes:
	    aa (list): Description
	    acc (list): Description
	    alpha (list): Description
	    antipar (list): Description
	    bp1 (list): Description
	    bp2 (list): Description
	    h_nho1 (list): Description
	    h_nho2 (list): Description
	    h_ohn1 (list): Description
	    h_ohn2 (list): Description
	    kappa (list): Description
	    moltyp (list): Description
	    num (list): Description
	    phi (list): Description
	    psi (list): Description
	    resnum (list): Description
	    sheet (list): Description
	    struct (list): Description
	    tco (list): Description
	    xyzca (list): Description
	"""
	
	def __init__(self):
			"""Summary
			"""
			self.num		= []
			self.resnum = []
			self.moltyp = []
			self.aa	  = []
			self.struct = []
			self.bp1		= []
			self.bp2		= []
			self.sheet		= []
			self.acc		= []
			self.h_nho1 = []
			self.h_ohn1 = []
			self.h_nho2 = []
			self.h_ohn2 = []
			self.tco		= []
			self.kappa  = []
			self.alpha  = []
			self.phi		= []
			self.psi		= []
			self.xyzca		= []
			self.antipar = []
			

	def parseDSSP(self, input_handle):
		"""Summary
		
		Args:
		    input_handle (TYPE): Description
		"""
		#input_handle = open(file, 'r')

		line_num = 0
		start=False
		for line in input_handle:
	
			if( re.search('#', line) ):
				start=True
				continue

			if( start and line[12:14].strip()!="!"):
				self.num.append(		line[0:5].strip() )
				self.resnum.append( line[5:10].strip() )
				self.moltyp.append( line[10:12].strip() )
				self.aa.append(	  line[12:14].strip() )
				self.struct.append( line[16:17].strip() )
				self.antipar.append( line[23:25] )
				self.bp1.append(		line[25:29].strip() )
				self.bp2.append(		line[29:33].strip() )
				self.sheet.append(		line[33].strip() )
				self.acc.append(		line[34:38].strip() )
				self.h_nho1.append( line[38:50].strip() )
				self.h_ohn1.append( line[50:61].strip() )
				self.h_nho2.append( line[61:72].strip() )
				self.h_ohn2.append( line[72:83].strip() )
				self.tco.append(		line[83:91].strip() )
				self.kappa.append(  line[91:97].strip() )
				self.alpha.append(  line[97:103].strip() )
				self.phi.append(		line[103:109].strip() )
				self.psi.append(		line[109:115].strip() )
				self.xyzca.append(		[float(line[115:122].strip()),float(line[122:129].strip()),float(line[129:136].strip())] )
				
		#input_handle.close()
	def getNums(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.num
	def getResnums(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.resnum
	def getChainType(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.moltyp
	def getAAs(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.aa
	def getSecStruc(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.struct
	def getAntiPar(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.antipar
	def getBP1(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.bp1
	def getBP2(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.bp2
	def getSheet(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.sheet
	def getACC(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.acc
	def getH_NHO1(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.h_nho1
	def getH_NHO2(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.h_nho2
	def getH_OHN1(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.h_ohn1
	def getH_OHN2(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.h_ohn2
	def getTCO(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.tco
	def getKAPPA(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.kappa
	def getALPHA(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.alpha
	def getPHI(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.phi
	def getPSI(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.psi
	def getXYZ(self):
		"""Summary
		
		Returns:
		    TYPE: Description
		"""
		return self.xyzca

def refinementCOMBO(cd, scores):
	"""Summary
	
	Args:
	    cd (TYPE): Description
	    scores (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	# cd is a matrix with nine spectra from each method.
	flag=0

	while(flag==0):
		cd = np.asarray(cd)
		scores=scores
		meancd= np.mean(cd, axis=0)
		stdcd=np.std(cd, axis=0)
		newcd=[]
		newscores=[]
		newq=[]
		qlist=list(range(0, len(cd)))
		for x,s,q in zip(cd, scores, qlist):
			bad=0
			for idx,y in enumerate(x):
				if(np.abs(y)>(np.abs(meancd[idx])+(1.5*stdcd[idx]))):
					bad+=1
			if(bad/len(cd[0]) <= (1./6.)):
				newcd.append(x)
				newscores.append(s)
				newq.append(q)
		if(len(newcd)==len(cd)):
			cd=np.asarray(newcd)
			scores=newscores
			qlist=newq
			flag=1
		else:
			cd=np.asarray(newcd)
			scores=newscores
			qlist=newq

	return(cd, scores, qlist)

def euclidDist(x, y):
	"""Summary
	
	Args:
	    x (TYPE): Description
	    y (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	distance=sqrt(sum((x-y)**2))
	return(distance)
   
def unit_vector(vector):
	"""Returns the unit vector of the vector.  
	
	Args:
	    vector (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
	"""Returns the angle in degrees between vectors 'v1' and 'v2'::
	
	>>> angle_between((1, 0, 0), (0, 1, 0))
	1.5707963267948966
	>>> angle_between((1, 0, 0), (1, 0, 0))
	0.0
	>>> angle_between((1, 0, 0), (-1, 0, 0))
	3.141592653589793
	
	Args:
	    v1 (TYPE): Description
	    v2 (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	return math.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))

def get_vector(p1, p2):
	"""Returns vector given two 3D coordinates 
	
	Args:
	    p1 (TYPE): Description
	    p2 (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	v=p1-p2

	return v

def distance_points(p1, v1, p2, v2):
	"""Returns close to the minimum distance between 2 lines that define secondary structure elements
	by using parametric rep of line and iterating t between 0 and 1 incrementing by 0.01.
	create array for each line with points alone line. Find distance between points. Take minimum 
	
	Args:
	    p1 (TYPE): Description
	    v1 (TYPE): Description
	    p2 (TYPE): Description
	    v2 (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	l1=[]
	l2=[]
	for i in range(0, 100):
		t=i*0.01

		# make array of l1 points
		l1.append(p1-v1*t)
		# make array of l2 points
		l2.append(p2-v2*t)

	l1=np.asarray(l1)
	l2=np.asarray(l2)

	mindistance=np.min(cdist(l1,l2, 'euclidean'))

	return(mindistance)

def distance_points_max(p1, v1, p2, v2):
	"""Summary
	
	Args:
	    p1 (TYPE): Description
	    v1 (TYPE): Description
	    p2 (TYPE): Description
	    v2 (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	l1=[]
	l2=[]
	for i in range(0, 100):
		t=i*0.01

		# make array of l1 points
		l1.append(p1-v1*t)
		# make array of l2 points
		l2.append(p2-v2*t)

	l1=np.asarray(l1)
	l2=np.asarray(l2)

	mindistance=np.max(cdist(l1,l2, 'euclidean'))

	return(mindistance)

def svd_line(data):
	"""Summary
	
	Args:
	    data (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	datamean = data.mean(axis=0)
	uu, dd, vv = np.linalg.svd(data - datamean)
	linepts = vv[0] * np.mgrid[-7:7:2j][:, np.newaxis]
	linepts += datamean

	return(linepts)

def rmsd(exp, pred):
	"""Summary
	
	Args:
	    exp (TYPE): Description
	    pred (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	total=0
	for i in range (0, len(pred)):
		x=(pred[i]-exp[i])**2
		total+=x
	result=sqrt(total/len(pred))

	return(result)

def getSS_multi(pdbfiles, directory, twistcutoff, mode, an):
	"""Summary
	
	Args:
	    pdbfiles (TYPE): Description
	    directory (TYPE): Description
	    twistcutoff (TYPE): Description
	    mode (TYPE): Description
	    an (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	all_ss=[]
	lstsq_ss=[]
	dssp_files=[]
	failed=[]
	all_HEO_full={}

	missing_dict={}
	# need to add this in

	with open(os.path.join(directory, "missing.txt"), "w") as missing_out:
		missing_out.write("")


	for i in range(0, len(pdbfiles)):
		
		try:
			
			f=pdbfiles[i]

			missing=0
			
			if(f in missing_dict):
				missing+=missing_dict[f]
				""" mmdict=MMCIF2Dict(os.path.join(outdir, pc+".cif"))
				if("_pdbx_unobs_or_zero_occ_residues.auth_seq_id" in mmdict):
					missing+=len(mmdict["_pdbx_unobs_or_zero_occ_residues.auth_seq_id"]) """
			else:
				with open(os.path.join(directory, f)) as fin:
					for l in fin.readlines():
						l=l.strip()
						missing_bool=re.search("^REMARK\s465\s+[A-Z]{3}\s[A-Z]\s*", l)
						if(missing_bool):
							missing+=1
			with open(os.path.join(directory, "missing.txt"), "a") as missing_out:
				missing_out.write(str(f) +" "+str(missing)+"\n")

			if(mode=="file"):
				dssp_input,success=pdbfile_to_dssp_local(directory,f)
				dssp_files.append(f)
			
			if(success==True):
				dd_ob=DSSPData()
				dd_ob.parseDSSP(dssp_input)
				cur_ch=99999
				heo_full=[]
				heo_ch=[]
				A=0 #right
				P=0
				H=0
				I=0
				G=0
				B=0
				T=0
				O=missing
				l=len(dd_ob.getSecStruc())+missing
				if(l==0):
					sys.exit()
				for ss, aorp, phi,ch in zip(dd_ob.getSecStruc(), dd_ob.getAntiPar(), dd_ob.getPHI(), dd_ob.getChainType()):
					if(ch!=cur_ch):
						cur_ch=ch
						if(len(heo_ch)>0):
							heo_full.append(heo_ch)
						heo_ch=[]

					if(ss=="H"):
						H+=1
						heo_ch.append("H")
					elif(ss=="E"):
						if(aorp[0].isupper() or aorp[1].isupper()):
							A+=1
						else:
							P+=1
						heo_ch.append("E")
					elif(ss=="I"):
						I+=1
						heo_ch.append("C")
					elif(ss=="G"):
						G+=1
						heo_ch.append("H")
					elif(ss=="B"):
						B+=1
						heo_ch.append("C")
					elif(ss=="T"):
						T+=1
						heo_ch.append("C")
					else:
						O+=1
						heo_ch.append("C")

				perp_count=getETOP_hbond(dd_ob, an)

				# for areas, name in zip(areas_all, name_all):
				# 	if(name!=False):
				# 		a_dists.append([areas])
				# 		names.append(name)
				heo_full.append(heo_ch)
				all_HEO_full[f]=(heo_full)
				A2=A*perp_count
				A1=A-A2
				
				all_ss.append([(H/l),G/l,A1/l,A2/l,P/l,T/l,I/l,B/l,O/l])
				lstsq_ss.append([(H/l),G/l,A1/l,A2/l,P/l,T/l,I/l+B/l+O/l])

				percentage=str(50+int(i/len(pdbfiles)*50))

			else:
				failed.append(pdbfiles[i])
				percentage=str(50+int(i/len(pdbfiles)*50))

				
		except:
			failed.append(pdbfiles[i])
			percentage=str(50+int(i/len(pdbfiles)*50))

			
		
	return(all_ss, lstsq_ss, dssp_files, failed, all_HEO_full)


def getETOP_hbond(dd_ob, an):
	"""
	curves following h-bonding patterns of beta residues
	
	Args:
	    dd_ob (TYPE): Description
	    an (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	ss=np.asarray(dd_ob.getSecStruc())
	xyz=np.asarray(dd_ob.getXYZ())
	res=np.asarray(dd_ob.getNums())
	ap=np.asarray(dd_ob.getAntiPar())
	bp1=np.asarray(dd_ob.getBP1())
	bp2=np.asarray(dd_ob.getBP2())
	sheet=np.asarray(dd_ob.getSheet())

	insheet="INIT"
	strands={}
	strand=[]
	strand_counter=0
	total_APE=0

	# make dictionary of resnum to xyz

	num_xyz_dict = dict(zip(res, xyz))


	# graph rep of h-bonded backbone vectors
	h_g = nx.Graph()
	

	for i in range(0, len(ss)):
		if(ss[i]=="E"):
			
			if(ap[i][0].isupper() or ap[i][1].isupper()):
				if(i<len(ss)-1):
					if(ss[i+1]=="E"):

						a1=num_xyz_dict[res[i]][0]-num_xyz_dict[res[i+1]][0]
						b1=num_xyz_dict[res[i]][1]-num_xyz_dict[res[i+1]][1]
						c1=num_xyz_dict[res[i]][2]-num_xyz_dict[res[i+1]][2]
						nodename=[a1,b1,c1]
						if(bp1[i] !="0"):
							a2=num_xyz_dict[bp1[i]][0]-num_xyz_dict[str(int(bp1[i])+1)][0]
							b2=num_xyz_dict[bp1[i]][1]-num_xyz_dict[str(int(bp1[i])+1)][1]
							c2=num_xyz_dict[bp1[i]][2]-num_xyz_dict[str(int(bp1[i])+1)][2]
							nodename2=[a2,b2,c2]
							h_g.add_node(json.dumps(nodename))
							h_g.add_node(json.dumps(nodename2))
							h_g.add_edge(json.dumps(nodename2),json.dumps(nodename))
						
						if(bp2[i] !="0"):
							a2=num_xyz_dict[bp2[i]][0]-num_xyz_dict[str(int(bp2[i])-1)][0]
							b2=num_xyz_dict[bp2[i]][1]-num_xyz_dict[str(int(bp2[i])-1)][1]
							c2=num_xyz_dict[bp2[i]][2]-num_xyz_dict[str(int(bp2[i])-1)][2]
							nodename2=[a2,b2,c2]
							h_g.add_node(json.dumps(nodename))
							h_g.add_node(json.dumps(nodename2))
							h_g.add_edge(json.dumps(nodename2),json.dumps(nodename))
						
						
	edge_list=list(h_g.edges())				
								
	A2=0
	
	for e in edge_list:
		e1=json.loads(e[0])
		e2=json.loads(e[1])
		a1=e1[0]
		b1=e1[1]
		c1=e1[2]
		a2=e2[0]
		b2=e2[1]
		c2=e2[2]
		angle=angleTwoLines(a1,b1,c1,a2,b2,c2)
		angle=math.degrees(angle)
		if((angle > int(an[1]) or angle < int(an[0]))):
			A2+=1



	if(len(edge_list)>0):
		A2=A2/len(edge_list)
	else:
		A2=0



	return(A2)

def isNMR(file, directory):
	"""Summary
	
	Args:
	    file (TYPE): Description
	    directory (TYPE): Description
	
	Returns:
	    TYPE: Description
	"""
	if(file.endswith(".pdb") or file.endswith(".ent") ):
		p = PDBParser(QUIET=True)
		structure = p.get_structure("X", os.path.join(directory, file))
		if(len(structure)>1):
			nmr=True
			new_pdb_files=[]
			for i, model in enumerate(structure):
				name=file[:-4]+"_"+str(i)+".pdb"
				io=PDBIO()
				io.set_structure(model)
				io.save(os.path.join(directory, name))
				new_pdb_files.append(name)
			#pdb_files.remove(file)
			#pdb_files=pdb_files+new_pdb_files
		else:
			new_pdb_files=[file]
	elif(file.endswith(".cif")):
		p = MMCIFParser(QUIET=True)
		structure = p.get_structure("X", os.path.join(directory, file))
		if(len(structure)>1):
			nmr=True
			new_pdb_files=[]
			for i, model in enumerate(structure):
				name=file[:-4]+"_"+str(i)+".pdb"
				io=PDBIO()
				io.set_structure(model)
				io.save(os.path.join(directory, name))
				new_pdb_files.append(name)
				if(file in missing_dict):
					missing_dict[name]=missing_dict[file]
			#pdb_files.remove(file)
			#pdb_files=pdb_files+new_pdb_files
		else:
			new_pdb_files=[file]

	return(new_pdb_files)

"""
WORK STARTS HERE 

"""

"""
get all ref lstsq_ss and lincomb_ss from reference files

"""

dsspfile_refs=[]
pcdfiles_refs=[]

with open(os.path.join(os.getcwd(), "ref_names.txt"), "r") as namesin:
	inputlines=namesin.readlines()
	for line in inputlines:
		f=re.split(r',', line)
		dsspfile_refs.append(f[0].strip())
		pcdfiles_refs.append(f[1].strip())

ref_lstsq_ss=[]

with open(os.path.join(os.getcwd(), "ref_lstsq_ss.txt"), "r") as lstsqin:
	inputlines=lstsqin.readlines()
	for line in inputlines:
		f=re.split(r',', line)
		ref_lstsq_ss.append([float(e) for e in f])

ref_lincomb_ss=[]

with open(os.path.join(os.getcwd(), "ref_lincomb_ss.txt"), "r") as lincombin:
	inputlines=lincombin.readlines()
	for line in inputlines:
		f=re.split(r',', line)
		ref_lincomb_ss.append([float(e) for e in f])

ref_lstsq_ss=np.asarray(ref_lstsq_ss)
ref_lincomb_ss=np.asarray(ref_lincomb_ss)

"""
get ref set CD

"""

ref_cd=[]

with open(os.path.join(os.getcwd(), "ref_cd.txt"), "r") as cdin:
	inputlines=cdin.readlines()
	for line in inputlines:
		f=re.split(r',', line)
		ref_cd.append([float(e) for e in f])

ref_cd=np.asarray(ref_cd)



"""
Basis spectra for LLSqs

"""

ref_basis=[]

with open(os.path.join(os.getcwd(), "ref_basis.txt"), "r") as basisin:
	inputlines=basisin.readlines()
	for line in inputlines:
		f=re.split(r',', line)
		ref_basis.append([float(e) for e in f])


ref_basis=np.asarray(ref_basis)


# functions for different modes

def makePred(f, lc_ss, lq_ss):
	"""Summary
	
	Args:
	    f (TYPE): Description
	    lc_ss (TYPE): Description
	    lq_ss (TYPE): Description
	"""
	# lincomb first

	candidates=[]

	for idx, ref_ss in enumerate(ref_lincomb_ss):
		candidates.append(euclidean(ref_ss, lc_ss))

	candidates_sorted, candidate_idxs_sorted, ref_lincomb_ss_sorted = zip(*[(x, y,z) for x, y,z in sorted(zip(candidates, list(range(0, len(candidates))), ref_lincomb_ss))])

	# get refs closest in ss to the query - caveat up to half more alpha/beta, vice versa. if < half greater then remainder is less.

	ss_greater=[]
	idx_greater=[]
	ss_lesser=[]
	idx_lesser=[]

	# get largest secondary structure proportion out of H,a1,a2,a3,P

	main_ss=np.argmax(np.asarray(lc_ss)[0:-1])

	for ss_sort, idxs_sorted in zip(ref_lincomb_ss_sorted, candidate_idxs_sorted):
		if(lc_ss[main_ss]< ss_sort[main_ss]):
			ss_greater.append(ss_sort)
			idx_greater.append(idxs_sorted)
		else:
			ss_lesser.append(ss_sort)
			idx_lesser.append(idxs_sorted)

	if(len(ss_greater)>no_cand/2):
		greater_no=int(no_cand/2)
	else:
		greater_no=int(len(ss_greater))
	if(len(idx_lesser) >= int(no_cand-greater_no)):
		lesser_no=int(no_cand-greater_no)
	else:
		lesser_no=len(idx_lesser)
		greater_no=int(no_cand-len(idx_lesser))
	

	if(greater_no ==0):
		candidate_ss=ref_lincomb_ss[np.asarray(idx_lesser[0:lesser_no])]
		candidate_cd=ref_cd[np.asarray(idx_lesser[0:lesser_no])]
	elif(lesser_no ==0):
		candidate_ss=ref_lincomb_ss[np.asarray(idx_greater[0:greater_no])]
		candidate_cd=ref_cd[np.asarray(idx_greater[0:greater_no])]

	else:
		candidate_ss=np.asarray(ref_lincomb_ss[np.asarray(idx_greater[0:greater_no])].tolist()+ref_lincomb_ss[np.asarray(idx_lesser[0:lesser_no])].tolist())
		candidate_cd=np.asarray(ref_cd[np.asarray(idx_greater[0:greater_no])].tolist()+ref_cd[np.asarray(idx_lesser[0:lesser_no])].tolist())

	# now need to minimise obj func to obtain proportions of each ref to get spectra
	objFun= lambda p: (euclidean(np.asarray(lc_ss), (np.sum(np.asarray([p]*len(candidate_ss.T)).T*candidate_ss, axis=0))))
	bnds = list(repeat([0.0, 1.0],len(candidate_ss)))
	x0 = list(repeat(0.1,len(candidate_ss)))
	X=minimize(objFun, x0, bounds=bnds)

	prediction=[]

	for x, cand_cd in zip(X.x, candidate_cd):
		
		prediction.append(cand_cd*x)

	lc_pred=np.sum(np.asarray(prediction), axis=0)

	# Lstsq second

	lq_pred=np.sum(ref_basis*np.asarray(lq_ss), axis=1)

	# combo pred third

	cd_refined_pred2, scores_refined_pred2, q_refined_pred2=refinementCOMBO(np.asarray([lq_pred, lc_pred]), list(range(0, 2)))
	if(len(cd_refined_pred2)>1):
		combo_pred=np.mean(np.asarray(cd_refined_pred2), axis=0)
	else:
		combo_pred=cd_refined_pred2[0]

	return(combo_pred)


def writeOutput(pred, f, outdir, success):
	if(success):
		output_f = os.path.join(outdir, str(f)+"_cd.txt")
		curr_time = str(datetime.datetime.now().strftime("%y-%m-%d, %H:%M:%S"))
		with open(output_f, "w") as fout:
			fout.write("# PDBMD2CD Standalone\n")
			fout.write("# E.D. Drew and R.W. Janes\n")
			fout.write("# https://doi.org/10.1093/nar/gkaa296\n")
			fout.write("# "+curr_time+"\n")
			fout.write("# Wavelength (nm)\tPredicted CD signal (delta epsilon)\n")
			for i, p in enumerate(pred):
				fout.write(str(260-i)+"\t"+str(p)+"\n")
	else:
		output_f = os.path.join(outdir, str(f)+"_cd.txt")
		curr_time = str(datetime.datetime.now().strftime("%y-%m-%d, %H:%M:%S"))
		with open(output_f, "w") as fout:
			fout.write("# PDBMD2CD Standalone\n")
			fout.write("# E.D. Drew and R.W. Janes\n")
			fout.write("# https://doi.org/10.1093/nar/gkaa296\n")
			fout.write("# "+curr_time+"\n")
			fout.write("# Wavelength (nm)\tPredicted CD signal (delta epsilon)\n")
			fout.write("Calculation Failed\n")
			

def writeLog(nFiles, nFails, logpath):
	curr_time=str(datetime.datetime.now().strftime("%y-%m-%d, %H:%M:%S"))
	with open(logpath, "a") as fout:
		fout.write(curr_time+":\t"+str(nFails)+"\t"+str(nFiles)+"\n")


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("indir", help="Path to directory containing .pdb files.")
	parser.add_argument("outdir", help="Path to directory for output files.")
	args = parser.parse_args()

	if(args.outdir):
		if(os.path.exists(args.outdir.strip())):
			outdir = args.outdir
		else:
			print("outdir - "+args.outdir+" - path doesn't exist")

	if(args.indir):
		if(os.path.exists(args.indir.strip())):
			indir = args.indir
		else:
			print("indir - "+args.indir+" - path doesn't exist")
	try:
		if(indir and outdir):
			print("Running DSSP....\n")
			all_files = []
			for f in os.listdir(indir):
				if(f.endswith(".pdb")):
					# check NMR/Models here
					files = isNMR(f, indir)
					all_files+=files

			all_lincomb_ss,all_lstsq_ss, dssp_files, failed, all_HEO_full = getSS_multi(all_files, indir, twistcutoff, "file", an)

			print("{0} failed calculations, {1} successful\n".format(str(len(failed)), str(len(dssp_files))))
			
			for f in failed:
				writeOutput([], f, outdir, False) # unecessary but just in case

			print("Making Predictions....\n")
			all_preds = []
			for f,  lc_ss, lstsq_ss in zip(dssp_files, all_lincomb_ss, all_lstsq_ss):
				if(f in failed):
					writeOutput([], f, outdir, False) # unecessary but just in case
				else:
					combo_pred = makePred(f, lc_ss, lstsq_ss)
					all_preds.append(combo_pred)

				writeOutput(combo_pred, f, outdir, True)

			mean_pred = np.mean(np.asarray(all_preds), axis=0)

			writeOutput(mean_pred, "average", outdir, True)
			
			writeLog(len(failed), len(combo_pred), config.LOGPATH)
			print("Log written: "+config.LOGPATH)
			print("Please remember to send log file to r.w.janes@qmul.ac.uk periodically so we can track usage. Thanks!")
			print("Done! Exiting")

		else:
			print("Exiting")
	except:
		print("Exiting")