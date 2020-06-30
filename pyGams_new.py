# Can be used if log file should be written instead of displayed
#save_stdout = sys.stdout
#fh = open("log.txt","w")
#sys.stdout = fh

"""
IMPORTANT: INSTRUCTIONS!!!

Changes from original:
0) edited: def readfromxls
1) new: def proc_read
2) edited: m.plants definition
3) new: def split.columns
4) new: def demsup_read
5) edited: m.demandall definition
6) edited: m.supimall definition
7) edited: def readfromgdx

Newly defined functions are commented.
Edited data is left but with a # in front.
If you want to check changes, just look for the new definitions of the aforementioned changes.
"""

from __future__ import division
import os
import sys
import pandas as pd
import numpy as np
import time
from gams import *
import pyGamstools as pt
import networkx as nx
import random
from IPython.core.debugger import Tracer
import Validation

ws = GamsWorkspace(working_directory = str(os.getcwd()),debug=2)
m = pd.DataFrame
read = pd.DataFrame

class modelData(object):

	#def __init__(self,fileMain,fileRen,fileDem,params,scenario):
	def __init__(self,fileMain,params,scenario):
			
		self.dbco = ws.add_database()
		self.dbprocess = ws.add_database()
		self.dbtransport = ws.add_database()
		self.dbstorage = ws.add_database()
		self.dbdsm = ws.add_database()
		self.dbsites = ws.add_database()
		self.dbtime = ws.add_database()
		self.dbflags = ws.add_database()
		self.dbPTDF = ws.add_database()
		self.dbDCDF = ws.add_database()
		self.dbPSDF = ws.add_database()
		
		print "added database"
		
	def readfromxls(self,fileMain,params,scenario):	
	#def readfromxls(self,fileMain,fileRen,fileDem,params,scenario):	
		
		global m
		# Create DataFrame from Excel file
		xls = pd.ExcelFile(fileMain)
		m.commodities = Validation.ValidateCommoditySheet(xls.parse('Commodity', convert_float=False))

		def proc_read(xls_file):
			"""
            reads process sheet from the input file, loading the appropriate evrys values.
    
			Args:
				- xls: main excel sheet
			
			Returns:
                evrys-compatible pandas DataFrame with process data
			
			Note:
			    supports just single input, and single output plus CO2.
				If fed with MIMO data, last value will be the considered one 
			"""
			import pandas as pd
			import numpy as np
			import math
			from datetime import datetime

            #imports Process-Commodity sheet, which contains the informations needed
			# for the creation of the evrys-compatible dataframe
			com = xls_file.parse('Process-Commodity', convert_float=False)
			expected_column_labels_2 = ['Process', 'Commodity', 'Direction', 'ratio', 'ratio-min']
    		Validation.CheckColumnNames(expected_column_labels_2, com, "Process-Commodity")
            #adds required evrys columns
			new_com_columns=['CoIn','CoOut','ratio_in','ratio_out','ratio_co2','ratio_in_min','ratio_out_min','eff','eff_min','cotwo']
			new_com = pd.DataFrame(columns=new_com_columns)
			new_com_rows=set(com['Process'])
            #adds string value to create dummy entries in the dataframe, which will be later 
			# changed for the appropriate values
			for i in range(len(new_com_columns)):
				this_column = new_com.columns[i]
				new_com[this_column] = ['n']*len(new_com_rows)
			new_com['Process']=new_com_rows
			new_com.set_index('Process', inplace=True)
            #adds the values which are already available from the excel sheet,
			# taking into consideration the direction of the process
			for i in range(len(com)):
				comprocess=com['Process'][i]
				if com['Direction'][i]=='In':
					new_com.at[comprocess,'CoIn']=com['Commodity'][i]
					#com.at[i,'CoOut']='Elec'
					new_com.at[comprocess,'ratio_in']=com['ratio'][i]
					new_com.at[comprocess,'ratio_in_min']=com['ratio-min'][i]
				elif com['Direction'][i]=='Out' and com['Commodity'][i]!='CO2':
					#com.at[i,'CoIn']='Elec'
					new_com.at[comprocess,'CoOut']=com['Commodity'][i]
					new_com.at[comprocess,'ratio_out']=com['ratio'][i]
					new_com.at[comprocess,'ratio_out_min']=com['ratio-min'][i]
				elif com['Direction'][i]=='Out' and com['Commodity'][i]=='CO2':
					new_com.at[comprocess, 'ratio_co2']=com['ratio'][i]
				else:
					print 'ValueError'
            #calculates the values which are needed for evrys from the values previously
			# loaded
			for row in list(new_com.index):
				new_com.at[row, 'eff']=new_com['ratio_out'][row]/new_com['ratio_in'][row]
				if new_com['ratio_co2'][row]=='n':
					new_com.at[row,'cotwo']=0
				else:
					new_com.at[row,'cotwo']=new_com['ratio_co2'][row]/new_com['ratio_in'][row]
				in_bool=math.isnan(new_com['ratio_in_min'][row])
				out_bool=math.isnan(new_com['ratio_out_min'][row])
				if in_bool==False and out_bool==False:
					new_com.at[row, 'eff_min']=new_com['ratio_out_min'][row]/new_com['ratio_in_min'][row]
				elif in_bool==True and out_bool==False:
					new_com.at[row, 'eff_min']=new_com['ratio_out_min'][row]/new_com['ratio_in'][row]
				elif in_bool==False and out_bool==True:
					new_com.at[row, 'eff_min']=new_com['ratio_out'][row]/new_com['ratio_in_min'][row]
				else:
					new_com.at[row, 'eff_min']=new_com['ratio_out'][row]/new_com['ratio_in'][row]
			new_com.drop(['ratio_in','ratio_out','ratio_in_min','ratio_out_min','ratio_co2'], axis=1, inplace=True)
			#loads the Process sheet, which is the one with the required evrys structure
			pro = xls_file.parse('Process', convert_float=False)
			#Check the column labels
    		expected_column_labels = ['Site', 'Process', 'inst-cap', 'cap-lo', 'cap-up', 'max-grad',
			'min-fraction', 'inv-cost', 'fix-cost', 'var-cost', 'wacc', 'y',
			'area-per-cap', 'act-up', 'on-off', 'start-cost', 'reserve-cost', 'ru',
			'rd', 'rumax', 'rdmax', 'cotwo', 'detail', 'lambda', 'heatmax',
			'maxdeltaT', 'heatupcost', 'su', 'sd', 'hotstart', 'pdt', 'pot',
			'prepow', 'pretemp', 'preheat', 'prestate', 'precaponline', 'year']
    		Validation.CheckColumnNames(expected_column_labels, pro ,'Process')

            #adds the required new columns to the Process dataframe
			new_columns=['CoIn','CoOut','eff','eff_min','cotwo']
			for column in new_columns:
				pro[column]=np.nan
				pro[column]=pro[column].fillna('n')
            #adds the correct values to the new dataframe, considering the fact that
			# in this new dataframe, process can appear multiple times, while in the 
			# previous data structure this didn't happen. Processes appeared just
			# one time each.
			for proc in new_com.index:
				for title in new_columns:
					n_appar=len(pro[title][proc])
					if n_appar==1:
						pro.at[proc, title]=new_com[title][proc]
					else:
						serie_ind=[proc]*n_appar
						serie_cont=[new_com[title][proc]]*n_appar
						serie=pd.Series(serie_cont, index=serie_ind)
						pro.at[proc, title]=serie
			pro.reset_index(inplace=True)
			pro.set_index(['Site','Process','CoIn','CoOut'], inplace=True)

			return pro

		#m.plants = xls.parse('Process', index_col=[0,1,2,3], convert_float=False)
		m.plants = proc_read(xls)
		
		m.transport = Validation.ValidateTransmissionSheet(xls.parse('Transmission', convert_float=False)) 
		m.storage   = Validation.ValidateStorageSheet(xls.parse('Storage', convert_float=False)) 
		m.dsm       = Validation.ValidateDsmSheet(xls.parse('DSM', convert_float=False)) 
		m.sites     = Validation.ValidateSiteSheet(xls.parse('Sites', convert_float=False))
		Validation.ValidateSiteNames(sites,commodities,plants,transport,storage,dsm) 

		def split_columns(columns, sep='.'):
			"""Split columns by separator into MultiIndex.

			Given a list of column labels containing a separator string (default: '.'),
			derive a MulitIndex that is split at the separator string.

			Args:
				- columns: list of column labels, containing the separator string
				- sep: the separator string (default: '.')

			Returns:
				a MultiIndex corresponding to input, with levels split at separator

			Example:
				>>> split_columns(['DE.Elec', 'MA.Elec', 'NO.Wind'])
				MultiIndex(levels=[['DE', 'MA', 'NO'], ['Elec', 'Wind']],
						labels=[[0, 1, 2], [0, 0, 1]])

			"""
			if len(columns) == 0:
				return columns
			column_tuples = [tuple(col.split('.')) for col in columns]
			return pd.MultiIndex.from_tuples(column_tuples)
        
        def demsup_read(request,xls_file):
			"""
            Reads the demand and supply data, and converts it from absolute value to 
			normalized value (sum of demand/supply for a certain site and commodity 
			is equal to 1.)

			Args:
				- request: which data is needed. This has to correspond to the sheet 
				  name in the fileMain excel sheet.
				- xls: main excel sheet
			
			Returns:
                evrys-compatible pandas DataFrame with demand/supply data
			"""
			demsup_temp = xls_file.parse(request).set_index(['t'])
			demsup=pd.DataFrame()
			#the next for loop is for the normalization of the values
			for index in list(demsup_temp.columns):
				demsup[index]=demsup_temp[index]/demsup_temp[index].sum()
			del demsup_temp
			#from here on the code is reading the demand input data
			# into a evrys-readable DataFrame
			demsup.columns = split_columns(demsup.columns, '.')
			demsup = demsup.stack().stack()
			demsup = demsup.to_frame()
			demsup = demsup.reset_index()
			#'level_1', 'level_2' and 0 are the column names which pandas automatically 
			# assigns to the DataFrame. Here the code is changing them to suitable 
			# column names.
			demsup = demsup.rename(columns = {'level_1':'co','level_2':'sit',0:'value'})
			#With the stack functions used before, the columns are not in order. That's
			# why the reindex function here below.
			column_ordered_list = ['t','sit','co','value']
			demsup = demsup.reindex(columns = column_ordered_list)
			demsup = demsup.set_index(['t','sit','co'])
			return demsup
		
		
       
		#m.demandall = read.from_csv('data/demand.csv', index_col=[0,1,2])
        m.demandall = demsup_read('Demand', fileMain)
		#m.supimall = read.from_csv('data/supim.csv', index_col=[0,1,2])
		m.supimall = demsup_read('SupIm', fileMain)	
		print "read"
		
		# Adapt data to scenario
		m = scenario(m)
		print "adapted to scenario"
		
		# Commodities
		co = self.dbco.add_set("co", 1, "commodities")
		for i in m.commodities.index.levels[1].unique()[:]:
			co.add_record(i.encode())
		co_sup_intermittent = self.dbco.add_set("co_sup_intermittent", 2, "intermittent commodities")
		co_sup_fix = self.dbco.add_set("co_sup_fix", 2, "fix constant commodities")
		co_sup_stock = self.dbco.add_set("co_sup_stock", 2, "stock commodities")
		co_dem = self.dbco.add_set("co_dem", 2, "demand commodities")
		for row in m.commodities.index.unique()[:]:
			myTuple=[str(x) for x in row]
			if m.commodities.loc[row]["type"]=="SupIm":
				co_sup_intermittent.add_record(myTuple)
			if m.commodities.loc[row]["type"]=="SUP-FIX":
				co_sup_fix.add_record(myTuple)
			if m.commodities.loc[row]["type"]=="Stock":
				co_sup_stock.add_record(myTuple)
			if m.commodities.loc[row]["type"]=="Demand":
				co_dem.add_record(myTuple)
		att_commodities = self.dbco.add_set("att_commodities", 1, "attributes of commodities")
		for i in m.commodities.columns.unique():
			att_commodities.add_record(i.encode())
		db_commodities = GamsParameter(self.dbco,"db_commodities", 3, "process attributes")
		for row in m.commodities.index:
			for att in m.commodities._get_numeric_data().columns:
				myTuple = [str(x) for x in row+tuple([att])]
				db_commodities.add_record(myTuple).value = m.commodities.loc[row][att.encode()]
		print "commodities"
		
		# Sites
		site = self.dbsites.add_set("site", 1, "sites")
		for i in m.sites.index.unique()[:]:
			site.add_record(i.encode())
		att_site = self.dbsites.add_set("att_sites", 1, "attributes of sites")
		for i in m.sites.columns.unique():
		    att_site.add_record(i.encode())
		db_site=GamsParameter(self.dbsites,"db_site",2,"site attributes")
		
		for row in m.sites.index:
			for att in m.sites._get_numeric_data().columns:
				myTuple1 = [str(x) for x in tuple([row])]
				myTuple2 = [str(x) for x in tuple([att])]
				myTuple = myTuple1+myTuple2
				db_site.add_record(myTuple).value = m.sites.loc[row][att.encode()]
		print "sites"
		
		# Processes
		pro=self.dbprocess.add_set("pro", 1, "plants")
		for i in m.plants.index.levels[1].unique()[:]:
			pro.add_record(i.encode())
		set_process = self.dbprocess.add_set("set_process", 4, "set process")
		for row in m.plants.index:
			myTuple = [str(x) for x in row]
			set_process.add_record(myTuple)	
		att_process = self.dbprocess.add_set("att_process", 1, "attributes")
		for k in m.plants.columns.unique():
			att_process.add_record(k.encode())
		db_process = GamsParameter(self.dbprocess, "db_process", 5, "process attributes")
		for row in m.plants.index:
			for att in m.plants.columns:
				myTuple = [str(x) for x in row+tuple([att])]
				db_process.add_record(myTuple).value = m.plants.loc[row][att.encode()]
		#import pdb; pdb.set_trace()
		print "process"
		
		# Storage
		sto = self.dbstorage.add_set("sto", 1, "storage")
		for i in m.storage.index.levels[1].unique()[:]:
			sto.add_record(i.encode())
		set_storage=self.dbstorage.add_set("set_storage", 3, "set storage")
		for row in m.storage.index:
			myTuple = [str(x) for x in row]
			set_storage.add_record(myTuple)	
		att_storage = self.dbstorage.add_set("att_storage",1,"attributes")
		for k in m.storage.columns.unique():
			att_storage.add_record(k.encode())
		db_storage = GamsParameter(self.dbstorage,"db_storage",4,"storage attributes")
		for row in m.storage.index:
			for att in m.storage.columns:
				myTuple = [str(x) for x in row+tuple([att])]
				db_storage.add_record(myTuple).value = m.storage.loc[row][att.encode()]
		print "storage"
		
		# DSM
		dsm = self.dbdsm.add_set("dsm", 1, "DSM")
		for i in m.dsm.index.levels[1].unique()[:]:
			dsm.add_record(i.encode())
		set_dsm=self.dbdsm.add_set("set_dsm", 2, "set dsm")
		for row in m.dsm.index:
			myTuple = [str(x) for x in row]
			set_dsm.add_record(myTuple)	
		att_dsm = self.dbdsm.add_set("att_dsm",1,"attributes")
		for k in m.dsm.columns.unique():
			att_dsm.add_record(k.encode())
		db_dsm = GamsParameter(self.dbdsm,"db_dsm",3,"dsm attributes")
		for row in m.dsm.index:
			for att in m.dsm.columns:
				myTuple = [str(x) for x in row+tuple([att])]
				db_dsm.add_record(myTuple).value = m.dsm.loc[row][att.encode()]
		dsm_past_up = GamsParameter(self.dbdsm,'dsm_past_up',3,'past upward dsm')
		for i in range(params['tStart']-int(max(m.dsm['delay'])), params['tStart']+params['nrSteps']+int(max(m.dsm['recovery']))-1):
			for row in m.dsm.index:
				myTuple = [str(i), str(row[0]), str(row[1])]
				dsm_past_up.add_record(myTuple).value = 0
		dsm_past_do = GamsParameter(self.dbdsm,'dsm_past_do',4,'past downward dsm')
		for i in range(params['tStart']-int(max(m.dsm['delay'])), params['tStart']+params['nrSteps']+int(max(m.dsm['recovery']))-1):
			for ii in range(i-3,i+4):
				if (ii>=params['tStart']-max(m.dsm['delay']) and ii<params['tStart']+params['nrSteps']+max(m.dsm['recovery'])-1):
					for row in m.dsm.index:
						myTuple = [str(i), str(ii), str(row[0]), str(row[1])]
						dsm_past_do.add_record(myTuple).value = 0
		print "DSM"
		
		# Timeseries
		t_ext = self.dbtime.add_set("t_ext", 1, "timesteps")	
		for i in range(params['tStart']-int(max(m.dsm['delay'])), params['tStart']+params['nrSteps']+int(max(m.dsm['recovery']))-1):
			t_ext.add_record(str(i))

		tm = self.dbtime.add_set("tm", 1, "timesteps")	
		for i in range(params['tStart'], params['tStart']+params['nrSteps']):
			tm.add_record(str(i))

		runcount = GamsParameter(self.dbtime, "runcount", 0, "runcount")
		runcount.add_record().value = params['runs']
		supim = GamsParameter(self.dbtime, "supim", 3, "supply intermittent")
		m.supim = m.supimall.loc[params['tStart']:params['tStart']+params['nrSteps']]
		m.supim = m.supim.reset_index()
		m.supim = m.supim.set_index(['t','sit','co'])
			
		# stochastische Alternative
		m.supim1=m.supim.copy()
		m.supim1.reset_index('t',inplace=True)
		for row in set(m.supim1.index):
			#Tracer()()
			deltam1=0
			for i in range(m.supim1.loc[row]["t"].min(),m.supim1.loc[row]["t"].max()+1):
				#Tracer()()
				row1=(i,)+row
				myTuple=[str(x) for x in row1]
				current=float(m.supim.loc[row1]["value"])
				if float(row1[0])>params['tStart']+params['nrSteps']-params['nrOverlap'] and params['stoch']==1 and row1[2]!="Hydro":
					#Tracer()()
					current1=random.normalvariate(current+deltam1, 0.02*((row1[0]-params['tStart'])**0.5)) 
					if current1>1:
						current1=1
					if current1<0:
						current1=0
					if current==0:
						current1=0
				else:
					current1=current
				#Tracer()()
				supim.add_record(myTuple).value=current1
				deltam1=current1-current
			
		demand = GamsParameter(self.dbtime, "demand", 3, "demand")
		m.demand = m.demandall.loc[params['tStart']:params['tStart']+params['nrSteps']]
		m.demand = m.demand.reset_index()
		m.demand = m.demand.set_index(['t','sit','co'])
		for row in m.demand.index:
			myTuple = [str(x) for x in row]
			demand.add_record(myTuple).value = float(m.demand.loc[row]["value"])
			#import pdb; pdb.set_trace()
		print "timeseries"
		
		# Transmission Grid
		m.transport.index.rename(['SitOut','SitIn','co','tr_type'], inplace=True)
		
		typ_transport = self.dbtransport.add_set("typ_transport", 1, "types of transmission lines")
		for k in m.transport.index.levels[3].unique():

			typ_transport.add_record(str(k))
				
		set_transport = self.dbtransport.add_set("set_transport", 4, "transport lines")
		for row in m.transport.index:
			myTuple=[str(x) for x in row]
			set_transport.add_record(myTuple)	
			
		att_transport = self.dbtransport.add_set("att_transport", 1, "attributes")
		for k in m.transport.columns.unique():
			att_transport.add_record(k.encode())
			
		db_transport=GamsParameter(self.dbtransport, "db_transport", 5, "transport attributes")
		#import pdb; pdb.set_trace()
		for row in m.transport.index:
			for att in m.transport.columns:
				myTuple=[str(x) for x in row+tuple([att])]
				db_transport.add_record(myTuple).value = m.transport.loc[row][att.encode()]
		print "grid"
		
		# # Time dependent startupcosts
		timestartcost=GamsParameter(self.dbtime,"timestartcost",5,"time dependent")


		for row in m.plants.index:
			myTuple=[str(x) for x in row]
			#import pdb; pdb.set_trace()
			tsc=np.zeros(params['nrSteps'])+m.plants.loc[tuple(myTuple)]["start-cost"]
			indexM=np.indices((1,params['nrSteps']))
			if (m.plants.loc[tuple(myTuple)]["heatupcost"]>0 and m.plants.loc[tuple(myTuple)]["act-lo"]>0  and m.plants.loc[tuple(myTuple)]["detail"]>3):
				tsc=tsc+m.plants.loc[tuple(myTuple)]["heatupcost"]*(1-np.exp(-m.plants.loc[tuple(myTuple)]["lambda"]*indexM[1,:]))
				#print '{0} at {1}'.format(tsc[0],myTuple)
                tsc=thinning(tsc[0][:],params['tol'])
                for i in range(1,params['nrSteps']):
					myTuple=[str(x) for x in row]
					timestartcost.add_record(tuple(["o"+str(i)])+tuple(myTuple)).value= tsc[i]
		
		print "startcost"
		
		# Flags
		set_flags=self.dbflags.add_set("set_flags",1,"plants")
		for i in m.flags.index.unique()[:]:
			set_flags.add_record(i.encode())
			
		db_flags=GamsParameter(self.dbflags,"db_flags",1,"flags")
		for i in m.flags.index:
			db_flags.add_record(str(i)).value=float(m.flags.loc[i]["value"])
		
		print "flags"
		
    #def readfromgdx(self,fileMain,fileRen,fileDem,params):
	def readfromgdx(self,fileMain,params):
		#Tracer()()
		self.dbflags = ws.add_database_from_gdx(os.path.join("%s/data" %str(os.getcwd()), "flags.gdx"))	
		self.dbtime=ws.add_database_from_gdx(os.path.join("%s/data" %str(os.getcwd()), "timeseries.gdx"))	
		self.dbco = ws.add_database_from_gdx(os.path.join("%s/data" %str(os.getcwd()), "commodity.gdx"))
		self.dbsites = ws.add_database_from_gdx(os.path.join("%s/data" %str(os.getcwd()), "sites.gdx"))			
		self.dbtransport = ws.add_database_from_gdx(os.path.join("%s/data" %str(os.getcwd()), "transport.gdx"))	
		self.dbstorage = ws.add_database_from_gdx(os.path.join("%s/data" %str(os.getcwd()), "storage.gdx"))	
		self.dbprocess = ws.add_database_from_gdx(os.path.join("%s/data" %str(os.getcwd()), "process.gdx"))
		#Tracer()()
		
		
	def set_flags(self,params):
		db_flags=self.dbflags.get_parameter('db_flags')
		db_flags.find_record("FTYPE").value=params['modelltype']
		db_flags.find_record("RMIP").value=params['relaxed']
		db_flags.find_record("TRANSTYPE").value=params['transtype']
		
		# Compileflags
		fobj_in = open("data/compileflags.dat",'w')
		fobj_in.write("$Setglobal FTYPE %s \n" %str(params['modelltype']))
		fobj_in.write("$Setglobal RMIP %s \n" %str(params['relaxed']))
		fobj_in.close()
		print "set flags"


	def write_data(self,params):
		#Tracer()()
		self.dbprocess.export(os.path.join("%s/data" %str(os.getcwd()), "process.gdx"))	
		self.dbtime.export(os.path.join("%s/data" %str(os.getcwd()), "timeseries.gdx"))	
		self.dbco.export(os.path.join("%s/data" %str(os.getcwd()), "commodity.gdx"))
		self.dbsites.export(os.path.join("%s/data" %str(os.getcwd()), "sites.gdx"))			
		self.dbtransport.export(os.path.join("%s/data" %str(os.getcwd()), "transport.gdx"))	
		self.dbstorage.export(os.path.join("%s/data" %str(os.getcwd()), "storage.gdx"))
		self.dbdsm.export(os.path.join("%s/data" %str(os.getcwd()), "dsm.gdx"))	

	
	def writeparams(self,params):
		self.dbflags.export(os.path.join("%s/data" %str(os.getcwd()), "flags.gdx"))	
		#Tracer()()
		if params['modelltype']<5:
			i=1
			fobj_in = open("cplex.opt",'w')
			fobj_in.write("preind 1\nlpmethod 4\nstartalg 4\nthreads 25\nmipemphasis 1\nheurfreq 5\nprobe 3\nreslim %s \noptcr %s" % (params['reslim'],params['mipgap']))
			fobj_in.close()

		if params['modelltype']==5:
			i=1
			fobj_in = open("cplex.opt",'w')
			fobj_in.write("preind 1\nlpmethod 4\nstartalg 4\nthreads 15\nusercutcall bchpython.inc mip=cplex")
			fobj_in.close()
		print "write data"
		
	def updatemodel(self,db_out,prestep,params):

		for pro in  self.dbprocess["set_process"]:
			prokeys=pro.keys
			self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"prepow"]).value=db_out["EPrOut"].find_record([str(prestep),prokeys[0],prokeys[1],prokeys[2],prokeys[3]]).level
			
			if self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"detail"]).value>1:
				self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"pretemp"]).value=db_out["Eth"].find_record([str(prestep),prokeys[0],prokeys[1],prokeys[2],prokeys[3]]).level
				self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"preheat"]).value=db_out["Eheat"].find_record([str(prestep),prokeys[0],prokeys[1],prokeys[2],prokeys[3]]).level
				
			if self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"detail"]).value>3:
				self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"prestate"]).value=db_out["StateOn"].find_record([str(prestep),prokeys[0],prokeys[1],prokeys[2],prokeys[3]]).level

			if self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"detail"]).value<4 and self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"detail"]).value>1:
				self.dbprocess["db_process"].find_record([prokeys[0],prokeys[1],prokeys[2],prokeys[3],"precaponline"]).value=db_out["CapOnl"].find_record([str(prestep),prokeys[0],prokeys[1],prokeys[2],prokeys[3]]).level	
		
		for sto in self.dbstorage["set_storage"]:
			stokeys=sto.keys
			self.dbstorage["db_storage"].find_record([stokeys[0],stokeys[1],stokeys[2],"precont"]).value=db_out["EstCon"].find_record([str(prestep),stokeys[0],stokeys[1],stokeys[2]]).level
			self.dbstorage["db_storage"].find_record([stokeys[0],stokeys[1],stokeys[2],"prepowout"]).value=db_out["EstOut"].find_record([str(prestep),stokeys[0],stokeys[1],stokeys[2]]).level
			self.dbstorage["db_storage"].find_record([stokeys[0],stokeys[1],stokeys[2],"prepowin"]).value=db_out["EstIn"].find_record([str(prestep),stokeys[0],stokeys[1],stokeys[2]]).level
			
		self.dbdsm['dsm_past_up'].clear()
		dsm_past_up = self.dbdsm.get_parameter('dsm_past_up')
		for dsm in db_out['dsm_up']:
			dsmkeys = dsm.keys
			dsm_past_up.add_record([dsmkeys[0],dsmkeys[1],dsmkeys[2]]).value = db_out['dsm_up'].find_record([dsmkeys[0],dsmkeys[1],dsmkeys[2]]).level
		
		self.dbdsm['dsm_past_do'].clear()
		dsm_past_do = self.dbdsm.get_parameter('dsm_past_do')
		for dsm in db_out['dsm_do']:
			dsmkeys = dsm.keys
			dsm_past_do.add_record([dsmkeys[0],dsmkeys[1],dsmkeys[2],dsmkeys[3]]).value = db_out['dsm_do'].find_record([dsmkeys[0],dsmkeys[1],dsmkeys[2],dsmkeys[3]]).level
		
		# Timeseries (have to be read from xls always)
		#Tracer()()
		self.dbtime['t_ext'].clear()
		t_ext = self.dbtime.get_set('t_ext')	
		for i in range(params['tStart']-int(max(m.dsm['delay'])), params['tStart']+params['nrSteps']+int(max(m.dsm['recovery']))-1):
			t_ext.add_record(str(i))

		self.dbtime["tm"].clear()
		tm=self.dbtime.get_set("tm")
		for i in range(params['tStart'],params['tStart']+params['nrSteps']):
			tm.add_record(str(i))

		self.dbtime["runcount"].clear()
		runcount=self.dbtime.get_parameter('runcount')
		runcount.add_record().value=params['runs']
		
		self.dbtime["supim"].clear()
		supim=self.dbtime.get_parameter('supim')
		m.supim=m.supimall.loc[params['tStart']:params['tStart']+params['nrSteps']]
		m.supim=m.supim.reset_index()
		m.supim=m.supim.set_index(['t','sit','co'])
		#for row in m.supim.index:
			#Tracer()()
		#	myTuple=[str(x) for x in row]
		#	supim.add_record(myTuple).value=float(m.supim.loc[row]["value"])
			
		# stochastische Alternative
		m.supim1=m.supim.copy()
		m.supim1.reset_index('t',inplace=True)
		for row in set(m.supim1.index):
			deltam1=0
			#Tracer()()
			for i in range(m.supim1.loc[row]["t"].min(),m.supim1.loc[row]["t"].max()+1):
				#Tracer()()
				row1=(i,)+row
				myTuple=[str(x) for x in row1]
				current=float(m.supim.loc[row1]["value"])
				if float(row1[0])>params['tStart']+params['nrSteps']-params['nrOverlap'] and params['stoch']==1 and row1[2]!="Hydro":
					#Tracer()()
					current1=random.normalvariate(current+deltam1, 0.02*((row1[0]-params['tStart'])**0.5)) 
					if current1>1:
						current1=1
					if current1<0:
						current1=0
					if current==0:
						current1=0
				else:
					current1=current
				supim.add_record(myTuple).value=current1
				deltam1=current1-current
			
		self.dbtime["demand"].clear()
		demand=self.dbtime.get_parameter('demand')
		m.demand=m.demandall.loc[params['tStart']:params['tStart']+params['nrSteps']]
		m.demand=m.demand.reset_index()
		m.demand=m.demand.set_index(['t','sit','co'])
		
		for row in m.demand.index:
			myTuple=[str(x) for x in row]
			demand.add_record(myTuple).value=float(m.demand.loc[row]["value"])
			
	def writePTDF(self,PTDF, nodelist, edgelist, firstrun):
	# Firstrun = 1 if writePTDF is called for the first time
	
		nodelistSlack=nodelist[:]
		#nodelistSlack.pop(0) commented as 0 should be written as well
		if firstrun==1:
			ItransAC=self.dbPTDF.add_set("ItransAC",3,"transport lines")
		self.dbPTDF["ItransAC"].clear()		
		ItransAC=self.dbPTDF.get_set("ItransAC")
		for row in edgelist:
			myTuple=[str(x) for x in row]
			myTuple=myTuple+['Elec']
			ItransAC.add_record(myTuple)	
		if firstrun==1:
			sits=self.dbPTDF.add_set("sits",1,"transport lines")	

		if firstrun==1:
			pPTDF=self.dbPTDF.add_parameter("pPTDF",4,"values")
		self.dbPTDF["pPTDF"].clear()	
		pPTDF=self.dbPTDF.get_parameter("pPTDF")
		e=0
		for row in edgelist:
			myTuple=[str(x) for x in row]
			myTuple=myTuple+['Elec']
			n=0
			for i in nodelistSlack:
				pPTDF.add_record(myTuple+[i.encode()]).value=PTDF[e,n]
				n=n+1
			e=e+1
		self.dbPTDF.export(os.path.join("%s/data" %str(os.getcwd()), "PTDF.gdx"))	

		return PTDF
		
	def writeDCDF(self,DCDF, edgelistAC, edgelistDC, firstrun):
	# Firstrun = 1 if writePTDF is called for the first time
	
		if firstrun==1:
			ItransDC=self.dbDCDF.add_set("ItransDC",3,"transport lines")
		self.dbDCDF["ItransDC"].clear()		
		ItransDC=self.dbDCDF.get_set("ItransDC")
		for row in edgelistDC:
			myTuple=[str(x) for x in row]
			myTuple=myTuple+['Elec']
			#print row
			ItransDC.add_record(myTuple)	

		if firstrun==1:
			pDCDF=self.dbDCDF.add_parameter("pDCDF",5,"values")
		self.dbDCDF["pDCDF"].clear()	
		pDCDF=self.dbDCDF.get_parameter("pDCDF")
		#Tracer()()
		e=0
		for row1 in edgelistAC:
			myTuple=[str(x) for x in row1]
			myTuple=myTuple+['Elec']
			n=0
			for row2 in edgelistDC:
				line2=[str(y) for y in row2]
				pDCDF.add_record(myTuple+line2).value=DCDF[e,n]
				n=n+1
			e=e+1
		self.dbDCDF.export(os.path.join("%s/data" %str(os.getcwd()), "DCDF.gdx"))	

		return DCDF
		
	def writePSDF(self,PSDF, edgelistAC, firstrun):
	# Firstrun = 1 if writePTDF is called for the first time
	
		if firstrun==1:
			pPSDF=self.dbPSDF.add_parameter("pPSDF",5,"values")
		self.dbPSDF["pPSDF"].clear()	
		pPSDF=self.dbPSDF.get_parameter("pPSDF")
		e=0
		for row1 in edgelistAC:
			myTuple=[str(x) for x in row1]
			myTuple=myTuple+['Elec']
			n=0
			for row2 in edgelistAC:
				line2=[str(y) for y in row2]
				pPSDF.add_record(myTuple+line2).value=PSDF[e,n]
				n=n+1
			e=e+1
		self.dbPSDF.export(os.path.join("%s/data" %str(os.getcwd()), "PSDF.gdx"))	

		return PSDF

def thinning(values,tolerance):
	if tolerance > 0:
		T = len(values)-1
		thinnedSteps = np.zeros((int(stepBound(values, tolerance))+1, 3))
		thinnedStepsMG=np.zeros(len(values))
		i = 1
		a = 1
		while a <= T:
			b = a
			while b+1 <= T and bestError(values[a], values[b+1]) < tolerance:
				b = b+1
			thinnedSteps[i,:] = [a,b,bestStep(values[a], values[b])]
			a = b + 1
			i = i+1
		if len(thinnedSteps)-1 >= i:
			thinnedSteps = thinnedSteps[:i]

			
		  #routine to get the thinned step in a format easier to handle for

		j=1
		for i in range(1,len(thinnedStepsMG)):
			if thinnedSteps[j,0]==i:
				thinnedStepsMG[i]=thinnedSteps[j,2]
				if j<len(thinnedSteps)-1:
					j=j+1
			else:
				thinnedStepsMG[i]=thinnedStepsMG[i-1]	
	else:
		thinnedStepsMG=values
	return thinnedStepsMG
          
	
def bestError(a,b):
	if a == 0 and b == 0:
		error = 0
	else:
		error = abs(a - b) / (a + b)
	return error
	
def bestStep(a,b):
	if a == 0 and b == 0:
		step = 0
	else:
		step = 2*a*b/(a + b);
	return step

def stepBound(values, tolerance):
	div = np.log( (1 + tolerance) / (1 - tolerance) )
	if values[1] == 0:
		if np.sum(values)==0:
			bound=1
		else:	
			bound = np.ceil(np.log(values[-1] / min(values[(values > 0)])) / div ) + 1
	else:
		bound = np.ceil(np.log(values[-1] / values[1] ) / div )
	return bound	

#fh.close()
# return to normal:
#sys.stdout = save_stdout


# get time

def write_output(var_list, marg_list, par_list, path_results, nrRuns):
	
	for output in var_list:
		for i in range(1,nrRuns+1):
			db_out = ws.add_database_from_gdx(path_results+'/ergebnis%s.gdx' % (i))
			#import pdb;
			#pdb.set_trace()
			out_var = pt.get_dataframe_var(db_out, output)
			if i==1:
				out_var_all = out_var


			else:

				#import pdb; pdb.set_trace()

				out_var_all = out_var.combine_first(out_var_all)
		out_var_all.to_csv(path_results+'/'+output.lower()+'all.csv')
	
	for output in marg_list:
		for i in range(1,nrRuns+1):
			db_out = ws.add_database_from_gdx(path_results+'/ergebnis%s.gdx' % (i))
			out_marg = pt.get_dataframe_margNodes(db_out, output)
			if i==1:
				out_marg_all = out_marg
			else:

				out_marg_all = out_marg.combine_first(out_marg_all)
		out_marg_all.to_csv(path_results+'/'+output.lower()+'all.csv')
		
	for output in par_list:
		for i in range(1,nrRuns+1):
			db_out = ws.add_database_from_gdx(path_results+'/ergebnis%s.gdx' % (i))
			out_par = pt.get_dataframe_par(db_out, output, i)
			if i==1:
				out_par_all = out_par
			else:
				out_par_all = out_par.combine_first(out_par_all)
		out_par_all.to_csv(path_results+'/'+output.lower()+'all.csv')
		
	
def makeGraph():
	GAC = nx.DiGraph()
	GDC = nx.DiGraph()

	# Nodes in order to define slacks
	Nodes = m.sites.copy()
	
	AC = m.transport.reset_index(level=[2, 3])
	AC = AC[AC.tr_type.map(lambda x: str(x).startswith('AC'))]
	AC = AC[['reactance']]
	AC.reset_index(inplace=True)	
	
	DC = m.transport.reset_index(level=[2, 3])
	DC = DC[DC.tr_type.map(lambda x: str(x).startswith('DC'))]
	
	DC = DC[['reactance']]
	DC.reset_index(inplace=True)	

	allnodes = pd.Series(Nodes.index.values)
	nodedictslack = Nodes.slacknode.to_dict()

	for y in allnodes.unique():
		GAC.add_node(y)
		GDC.add_node(y)

	for x in AC.values:
		GAC.add_weighted_edges_from([tuple(x)])

	for x in DC.values:
		GDC.add_weighted_edges_from([tuple(x)])

	nx.set_node_attributes(GAC,nodedictslack,'slack')
	nx.set_node_attributes(GDC,nodedictslack,'slack')

	return GAC, GDC
	
def sortKeyWith2ndThen1stListValue(item):
  return item[1], item[0]
	
def generatePTDF(G):

	nodelist= sorted(G.nodes())
	edgelist= list(G.edges())
	#edgelist.sort(key=sortKeyWith2ndThen1stListValue) # sortieren nach zweitem Node, anpassen an xls am besten
	A=nx.incidence_matrix(G,nodelist, edgelist,True)
	A=A.toarray()
	A=A.transpose() #andere Defintion in networkx
	A=np.asmatrix(A)
	A=-A

	X=np.zeros((len(edgelist),len(edgelist)))
	for i in range(0,len(edgelist)):
		X[i,i]=G.edges[edgelist[i][0], edgelist[i][1]]['weight']
	#import pdb;pdb.set_trace()
	slacks=nx.get_node_attributes(G,'slack')
	zipped=sorted(slacks.items(),key=getKey)
	unzipped=zip(*zipped)
	uzarray=np.asarray(unzipped[1])
	referencenodes=np.where(uzarray==1)[0]	
	B=X.copy()
	B[B!=0]=-1/B[B!=0] #invertiere Reaktancen to get susceptances
	B=np.asmatrix(B)
	Aux1=(B*A)
	Aux2=A.transpose()*B*A
	#delete slack node
	Aux1=np.delete(Aux1,referencenodes,axis=1) # column
	Aux2=np.delete(Aux2,referencenodes,axis=0) # row
	Aux2=np.delete(Aux2,referencenodes,axis=1) # column
	#import pdb;	pdb.set_trace()
	PTDF=Aux1*(np.linalg.inv(Aux2))
	refinsert=referencenodes.copy()
	for i in range(0,len(refinsert)):
		refinsert[i]=refinsert[i]-i
	PTDF=np.insert(PTDF,refinsert,0,axis=1)

	#import pdb;pdb.set_trace()
	return A, B, PTDF, nodelist, edgelist
	
	
def generateDCDF(G, PTDF):
	nodelist= sorted(G.nodes())
	edgelist= list(G.edges())
	edgelist.sort(key=sortKeyWith2ndThen1stListValue) # sortieren nach zweitem Node, anpassen an xls am besten
	A=nx.incidence_matrix(G,nodelist, edgelist,True)
	A=A.toarray()
	A=A.transpose() #andere Defintion in networkx
	A=np.asmatrix(A)
	A=-A
	DCDF=-PTDF*A.transpose()
	#import pdb; pdb.set_trace()
	#DCDF=np.insert(DCDF,0,0,axis=1)

	return DCDF, edgelist
	
	
def generatePSDF(G, PTDF):
	nodelist= sorted(G.nodes())
	edgelist= list(G.edges())
	edgelist.sort(key=sortKeyWith2ndThen1stListValue) # sortieren nach zweitem Node, anpassen an xls am besten
	A=nx.incidence_matrix(G,nodelist, edgelist,True)
	A=A.toarray()
	A=A.transpose() #andere Defintion in networkx
	A=np.asmatrix(A)
	A=-A
	# Read susceptances
	X=np.zeros((len(edgelist),len(edgelist)))
	#B=np.zeros((len(edgelist),len(edgelist)))
	for i in range(0,len(edgelist)):
		X[i,i]=G.edges[edgelist[i][0], edgelist[i][1]]['weight']

	B=X.copy()
	B[B!=0]=-1/B[B!=0] #invertiere Reaktancen to get susceptances
	B=np.asmatrix(B)
	
	Aux1=(B*A)
	Aux2=B

	PSDF=Aux2-PTDF*Aux1.transpose()
	#import pdb;pdb.set_trace()
	return PSDF


def getKey(item):
	return item[0]