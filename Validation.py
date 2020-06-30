import os
import pandas as pd
import future


def CheckColumnNames(expected_column_labels,dataframe,sheetname):
    """
    First Arg: expected labels of columns as a list
    Second Arg: dataframe of the sheet
    Third Arg: The name of the sheet
    """
    for label in expected_column_labels:
        if label not in dataframe.columns:
            raise KeyError("Column '{}' could not be found in the {} sheet. Please check the column name!".format(label,sheetname))    

def ValidateSiteSheet(dataframe):
    """Excel data is processed and checked for further use.
    """
    expected_column_labels=['Name', 'area', 'slacknode', 'lat', 'long', 'ctrarea', 'primpos',
       'primneg', 'secpos', 'secneg', 'terpos', 'terneg', 'syncharea',
       'htworegion']
    CheckColumnNames(expected_column_labels,dataframe,"Site")

    exclusive_urbs_columns = ["area"]
    dataframe = dataframe.drop(exclusive_urbs_columns, axis="columns")

    dataframe = dataframe.rename(columns={"Name":"Site"})
    dataframe = dataframe.set_index(["Site"])

    columns_ordered = ['slacknode', 'lat', 'long', 'ctrarea', 'primpos', 'primneg',
       'secpos', 'secneg', 'terpos', 'terneg', 'syncharea', 'htworegion']
    dataframe = dataframe.reindex(columns=columns_ordered)
    return dataframe

#Modifies the commodity sheet and relabels the data for further use.
def ValidateCommoditySheet(dataframe):
    """Excel data is processed and checked for further use.
    """
    expected_column_labels = ['Site', 'Commodity', 'Type', 'price', 'max', 'maxperhour', 'annual',
       'losses']
    CheckColumnNames(expected_column_labels,dataframe,'Commodity')

    exclusive_urbs_columns_commodity = ["max", "maxperhour"]
    dataframe = dataframe.drop(exclusive_urbs_columns_commodity, axis='columns')

    dataframe = dataframe.rename(columns={'Commodity':'Co',"Type":"type"})
    dataframe = dataframe.set_index(["Site","Co"])
    columns_ordered=["price","annual","losses","type"]
    dataframe = dataframe.reindex(columns=columns_ordered)
    return dataframe  

#Modifies the process sheet and relabels the data for further use. 
#Incomplete
def ValidateProcessSheet(dataframe):
    """Excel data is processed and checked for further use.
    """
    #Check the column labels
    expected_column_labels = ['Site', 'Process', 'inst-cap', 'cap-lo', 'cap-up', 'max-grad',
       'min-fraction', 'inv-cost', 'fix-cost', 'var-cost', 'wacc', 'y',
       'area-per-cap', 'act-up', 'on-off', 'start-cost', 'reserve-cost', 'ru',
       'rd', 'rumax', 'rdmax', 'cotwo', 'detail', 'lambda', 'heatmax',
       'maxdeltaT', 'heatupcost', 'su', 'sd', 'hotstart', 'pdt', 'pot',
       'prepow', 'pretemp', 'preheat', 'prestate', 'precaponline', 'year']
    CheckColumnNames(expected_column_labels, dataframe,'Process')

    expected_column_labels_2 = ['Process', 'Commodity', 'Direction', 'ratio', 'ratio-min']
    CheckColumnNames(expected_column_labels_2, dataframe  ,"Process-Commodity")
    
    # #Remove urbs exclusive columns
    # exclusive_urbs_columns_process = ["cap-lo","cap-up","max-grad","inv-cost","fix-cost","var-cost","wacc","y","area-per-cap"]
    # dataframe = dataframe.drop(exclusive_urbs_columns_process, axis='columns')

    #relabel columns
    dataframe = dataframe.rename(columns={"Process":"Pro", "min-fraction":"act-lo"})

    # #Set indexes 
    # dataframe = dataframe.set_index(["Site","Pro"]) #"CoIn","CoOut" is missing
    
    return dataframe

#Modifies the transmission sheet and relabels the data for further use. 
def ValidateTransmissionSheet(dataframe):
    """Excel data is processed and checked for further use.
    """
    #Check the columns
    expected_column_labels = ['Site In', 'Site Out', 'Transmission', 'Commodity', 'eff', 'inv-cost',
       'fix-cost', 'var-cost', 'inst-cap', 'cap-lo', 'cap-up', 'wacc',
       'depreciation', 'reactance', 'difflimit', 'base_voltage',
       'cap-up-therm', 'angle-up', 'u2b', 'dc-flow', 'length', 'react-pu',
       'PSTmax', 'idx']
    CheckColumnNames(expected_column_labels,dataframe,"Transmission")

    #Drop unnecessary columns
    exclusive_urbs_columns_transmission = ["inv-cost","fix-cost","wacc","depreciation","difflimit","base_voltage"]
    dataframe = dataframe.drop(exclusive_urbs_columns_transmission, axis="columns")

    #Relabel the columns
    dataframe = dataframe.rename(columns={"Site In":"SitIn", "Site Out":"SitOut", 
    "Commodity":"Co", "cap-lo":"act-lo", "cap-up":"act-up","Transmission":"tr_type"}) 

    #Set indexes
    dataframe = dataframe.set_index(["SitIn","SitOut","Co","tr_type"])

    #Ordering
    columns_ordered = ['SitIn', 'SitOut', 'Co', 'eff', 'var-cost', 'inst-cap', 'act-lo',
       'act-up', 'reactance', 'cap-up-therm', 'angle-up', 'u2b', 'dc-flow',
       'length', 'react-pu', 'tr_type', 'PSTmax', 'idx']
    dataframe = dataframe.reindex(columns=columns_ordered)

    return dataframe

def ValidateStorageSheet(dataframe):
    """Excel data is processed and checked for further use.
    """
    #Check Columns
    expected_column_labels = ['Site', 'Storage', 'Commodity', 'inst-cap-c', 'cap-lo-c', 'cap-up-c',
       'inst-cap-p', 'cap-lo-p', 'cap-up-p', 'eff-in', 'eff-out', 'inv-cost-p',
       'inv-cost-c', 'fix-cost-p', 'fix-cost-c', 'var-cost-p', 'var-cost-c',
       'wacc', 'depreciation', 'init', 'discharge', 'ep-ratio', 'inst-cap-pi',
       'inst-cap-po', 'var-cost-pi', 'var-cost-po', 'act-lo-pi', 'act-up-pi',
       'act-lo-po', 'act-up-po', 'act-lo-c', 'act-up-c', 'precont', 'prepowin',
       'prepowout', 'ru', 'rd', 'rumax', 'rdmax', 'seasonal', 'ctr']
    CheckColumnNames(expected_column_labels,dataframe,"Storage")

    #Drop unneccesary columns
    exclusive_urbs_columns = ["cap-lo-c","cap-up-c","inst-cap-p","cap-lo-p","cap-up-p","inv-cost-p","inv-cost-c","fix-cost-p",
    "fix-cost-c","var-cost-p","wacc","discharge","ep-ratio"]
    dataframe = dataframe.drop(exclusive_urbs_columns, axis='columns')

    #Relabel columns for further use
    dataframe = dataframe.rename(columns={"Storage":"Sto","Commodity":"Co"})

    #Indexing
    dataframe = dataframe.set_index(["Site","Sto","Co"])

    #Ordering
    column_names_order = ['inst-cap-pi', 'inst-cap-po', 'inst-cap-c',
       'eff-in', 'eff-out', 'var-cost-pi', 'var-cost-po', 'var-cost-c',
       'act-lo-pi', 'act-up-pi', 'act-lo-po', 'act-up-po', 'act-lo-c',
       'act-up-c', 'precont', 'prepowin', 'prepowout', 'ru', 'rd', 'rumax',
       'rdmax', 'seasonal', 'ctr']
    dataframe = dataframe.reindex(columns=column_names_order)

    return dataframe

def ValidateDsmSheet(dataframe):
    """Excel data is processed and checked for further use.
    """
    expected_column_labels = ['Site', 'Commodity', 'delay', 'eff', 'recov', 'cap-max-do',
       'cap-max-up', 'rel-inst-cap', 'var-cost']
    CheckColumnNames(expected_column_labels,dataframe,"DSM")

    #Drop unnecesary columns
    exclusive_urbs_columns = ["cap-max-do","cap-max-up"]
    dataframe = dataframe.drop(exclusive_urbs_columns, axis="columns")

    #Relabel columns for further use
    dataframe = dataframe.rename(columns={'Commodity':"Co","recov":"recovery"})
    
    #Indexing
    dataframe = dataframe.set_index(["Site","Co"])

    #ordering
    columns_ordered=['rel-inst-cap', 'eff', 'delay', 'recovery', 'var-cost']
    dataframe = dataframe.reindex(columns=columns_ordered)
    return dataframe

def ValidateSiteNames(dfSite,dfCommodity,dfProcess,dfTransmission,dfStorage,dfDSM):
    sites = dfSite.index.tolist()
    #Check site names in commodity sheet against site sheet
    for site in dfCommodity.index.levels[0].tolist():
        if site not in sites:
            raise KeyError("The site name '{}' in Commodity sheet is not listed in the sheet 'Site'!".format(site))
    #check site names in process sheet against site sheet
    for site in dfProcess.index.levels[0].tolist():
        if site not in sites:
            raise KeyError("The site name '{}' in Process sheet is not listed in the sheet 'Site'!".format(site))
    #Check site names in transmission sheet against site sheet
    for site in dfTransmission.index.levels[0].tolist():
        if site not in sites:
            raise KeyError("The site name '{}' in Transmission sheet at the column 'Site In' is not listed in the sheet 'Site'!".format(site))
    for site in dfTransmission.index.levels[1].tolist():
        if site not in sites:
            raise KeyError("The site name '{}' in Transmission sheet at the column 'Site Out' is not listed in the sheet 'Site'!".format(site))
    #Check site names in Storage sheet against site sheet
    for site in dfStorage.index.levels[0].tolist():
        if site not in sites:
            raise KeyError("The site name '{}' in Storage sheet is not listed in the sheet 'Site'!".format(site))
    #Check site names in DSM sheet against site sheet
    for site in dfDSM.index.levels[0].tolist():
        if site not in sites:
            raise KeyError("The site name '{}' in DSM sheet is not listed in the sheet 'Site'!".format(site))

# #Site
# excel="data\onesyncharea.xlsx"
# xls = pd.ExcelFile(excel)
# sites = xls.parse('Sites', index_col=[0], convert_float=False)
# print(sites.head())
# print("--^--"*42)

# newExcel="data\propens.xlsx"
# xls = pd.ExcelFile(newExcel)
# sites = xls.parse("Site", convert_float=False)
# sites = ValidateSiteSheet(sites)
# print(sites.head())

# #Commodity
# excel="data\onesyncharea.xlsx"
# xls = pd.ExcelFile(excel)
# commodities = xls.parse('Commodities', index_col=[0,1], convert_float=False)
# print(commodities.head())
# print('--^--'*42)

# newExcel="data\propens.xlsx"
# xls = pd.ExcelFile(newExcel)
# commodity = xls.parse('Commodity', convert_float=False)
# commodity = ValidateCommoditySheet(commodity)
# print(commodity.head())

#Process
excel = "data\onesyncharea.xlsx"
xls = pd.ExcelFile(excel)
plants = xls.parse('Process', index_col=[0,1,2,3], convert_float=False)
print(plants.head())
print('--^--'*42)

# newExcel="data\propens.xlsx"
# xls = pd.ExcelFile(newExcel)
# plants = xls.parse('Process', convert_float=False)
# plants = ValidateProcessSheet(plants)
# print(plants)

#Transmission
# excel="data\onesyncharea.xlsx"
# xls = pd.ExcelFile(excel)
# transport = xls.parse('Transmission', index_col=[0,1,2,15], convert_float=False)
# print(transport.head())
# print('--^--'*42)

# newExcel="data\propens.xlsx"
# xls = pd.ExcelFile(newExcel)
# transport = xls.parse('Pr', convert_float=False)
# transport = ValidateTransmissionSheet(transport)
# print(transport.head())

# #Storage
# excel="data\onesyncharea.xlsx"
# xls = pd.ExcelFile(excel)
# storage = xls.parse('Storage', index_col=[0,1,2], convert_float=False)
# print(storage.head())
# print('--^--'*42)

# newExcel="data\propens.xlsx"
# xls = pd.ExcelFile(newExcel)
# storage = xls.parse("Storage", convert_float=False)
# storage = ValidateStorageSheet(storage)
# print(storage.head())

# #DSM
# excel="data\onesyncharea.xlsx"
# xls = pd.ExcelFile(excel)
# dsm = xls.parse('DSM', index_col=[0,1], convert_float=False)
# print(dsm.head())
# print('--^--'*42)

# newExcel="data\propens.xlsx"
# xls = pd.ExcelFile(newExcel)
# dsm = xls.parse("DSM", convert_float=False)
# dsm = ValidateDsmSheet(dsm)
# print(dsm.head())


# newExcel="data\propens.xlsx"
# xls = pd.ExcelFile(newExcel)

# commodities = xls.parse('Commodity', convert_float=False)
# commodities = ValidateCommoditySheet(commodities)
# plants = xls.parse('Process', convert_float=False)
# plants = ValidateProcessSheet(plants)
# transport = xls.parse('Transmission', convert_float=False)
# transport = ValidateTransmissionSheet(transport)
# storage = xls.parse('Storage', convert_float=False)
# storage = ValidateStorageSheet(storage)
# dsm = xls.parse('DSM', convert_float=False)
# dsm = ValidateDsmSheet(dsm)
# sites = xls.parse('Site', convert_float=False)
# sites = ValidateSiteSheet(sites)

# # print(transport.index)

# ValidateSiteNames(sites,commodities,plants,transport,storage,dsm)


# newExcel="data\propens.xlsx"
# xls = pd.ExcelFile(newExcel)
# plants = xls.parse('Process-Commodity', convert_float=False)
# print(plants.head())  

# CoIn = {}
# plants = plants.set_index('Process')
# for plantname in plants['Process']:
#     print(plants.loc[plantname])
        

        




