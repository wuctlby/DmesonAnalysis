'''
python script with helper functions to load objects from task
'''

import sys
from ROOT import TFile  # pylint: disable=import-error,no-name-in-module

# pylint: disable=too-many-branches,too-many-statements, too-many-return-statements
def LoadSparseFromTask(infilename, inputCfg, no_List_isMC=False):
    '''
    Method to retrieve sparses from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file
    - flag to whether load 'isMC' and list from config
        - if True, load sparses for MC and no list
        - if False, load 'isMC' from config and load sparse from list
        - default is False

    Returns
    ----------
    - list of sparses with reconstructed quantities for all candidates,
      and prompt, FD D mesons and bkg candidates (only if MC)
    - list of sparses with generated quantities for prompt and FD D mesons (only if MC)
    '''
    print('Loading THnSparses from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None, None
    
    if not no_List_isMC:
        inlistData = indirData.Get(inputCfg['listname'])
        if not inlistData:
            print(f'List {inputCfg["listname"]} not found!')
            return None, None
        isMC = inputCfg['isMC'] if 'isMC' in inputCfg else False
    else:
        isMC = True

    sparses, sparsesGen = {}, {}
    if inputCfg['sparsenameAll']:
        if not no_List_isMC:
            sparses['RecoAll'] = inlistData.FindObject(inputCfg['sparsenameAll']) # not mandatory for MC
        else:
            sparses['RecoAll'] = indirData.Get(inputCfg['sparsenameAll'])
        if not sparses['RecoAll']:
            print(f'ERROR: sparse {inputCfg["sparsenameAll"]} not found!')
            return None, None
    if isMC:
        if ((inputCfg['sparsenamePrompt'] is not None and inputCfg['sparsenameAll'] == inputCfg['sparsenamePrompt']) or
                (inputCfg['sparsenameFD'] is not None and inputCfg['sparsenameAll'] == inputCfg['sparsenameFD'])):
            print('ERROR: do not use the same object for different sparses, this gives an error when merged! Exit')
            sys.exit()
        if inputCfg['sparsenamePrompt'] is not None:
            if not no_List_isMC:
                sparses['RecoPrompt'] = inlistData.FindObject(inputCfg['sparsenamePrompt'])
            else:
                sparses['RecoPrompt'] = indirData.Get(inputCfg['sparsenamePrompt'])
                sparses['RecoPrompt'].SetName("SparseRecoPrompt")
            if not sparses['RecoPrompt']:
                print(f'ERROR: sparse {inputCfg["sparsenamePrompt"]} not found!')
                return None, None
        if inputCfg['sparsenameFD'] is not None:
            if not no_List_isMC:
                sparses['RecoFD'] = inlistData.FindObject(inputCfg['sparsenameFD'])
            else:
                sparses['RecoFD'] = indirData.Get(inputCfg['sparsenameFD'])
                sparses['RecoFD'].SetName("SparseRecoFD")
            if not sparses['RecoFD']:
                print(f'ERROR: sparse {inputCfg["sparsenameFD"]} not found!')
                return None, None
        if not no_List_isMC:
            sparsesGen['GenPrompt'] = inlistData.FindObject(inputCfg['sparsenameGenPrompt'])
        else:
            sparsesGen['GenPrompt'] = indirData.Get(inputCfg['sparsenameGenPrompt'])
            sparsesGen['GenPrompt'].SetName("SparseGenPrompt")
        if not sparsesGen['GenPrompt']:
            print(f'ERROR: sparse {inputCfg["sparsenameGenPrompt"]} not found!')
            return None, None
        if not no_List_isMC:
            sparsesGen['GenFD'] = inlistData.FindObject(inputCfg['sparsenameGenFD'])
        else:
            sparsesGen['GenFD'] = indirData.Get(inputCfg['sparsenameGenFD'])
            sparsesGen['GenFD'].SetName("SparseGenFD")
        if not sparsesGen['GenFD']:
            print(f'ERROR: sparse {inputCfg["sparsenameGenFD"]} not found!')
            return None, None
        
        # For D0, the reflection and reconstruction are in the same sparse
        if inputCfg['enableRef']:
            if inputCfg['sparsenameRefl'] is not None:
                if not no_List_isMC:
                    sparses['RecoRefl'] = inlistData.FindObject(inputCfg['sparsenameRefl'])
                else:
                    sparses['RecoRefl'] = indirData.Get(inputCfg['sparsenameRefl'])
                    sparses['RecoRefl'].SetName("SparseRefl")
                    sparses['RecoReflPrompt'] = sparses['RecoRefl'].Clone('SparseReflPrompt')
                    sparses['RecoReflFD'] = sparses['RecoRefl'].Clone('SparseReflFD')
            if not sparses['RecoRefl']:
                print(f'ERROR: sparse {inputCfg["sparsenameRefl"]} not found!')
                return None, None

        if inputCfg['enableSecPeak']:
            if inputCfg['sparsenamePromptSecPeak'] is not None:
                if not no_List_isMC:
                    sparses['RecoSecPeakPrompt'] = inlistData.FindObject(inputCfg['sparsenamePromptSecPeak'])
                else:
                    sparses['RecoSecPeakPrompt'] = indirData.Get(inputCfg['sparsenamePromptSecPeak'])
                    sparses['RecoSecPeakPrompt'].SetName("SparseRecoSecPeakPrompt")
                if not sparses['RecoSecPeakPrompt']:
                    print(f'ERROR: sparse {inputCfg["sparsenamePromptSecPeak"]} not found!')
                    return None, None
            if inputCfg['sparsenameFDSecPeak'] is not None:
                if not no_List_isMC:
                    sparses['RecoSecPeakFD'] = inlistData.FindObject(inputCfg['sparsenameFDSecPeak'])
                else:
                    sparses['RecoSecPeakFD'] = indirData.Get(inputCfg['sparsenameFDSecPeak'])
                    sparses['RecoSecPeakFD'].SetName("SparseRecoSecPeakFD")
                if not sparses['RecoSecPeakFD']:
                    print(f'ERROR: sparse {inputCfg["sparsenameFDSecPeak"]} not found!')
                    return None, None
            if not no_List_isMC:
                sparsesGen['GenSecPeakPrompt'] = inlistData.FindObject(inputCfg['sparsenameGenPromptSecPeak'])
            else:
                sparsesGen['GenSecPeakPrompt'] = indirData.Get(inputCfg['sparsenameGenPromptSecPeak'])
                sparsesGen['GenSecPeakPrompt'].SetName("SparseGenSecPeakPrompt")
            if not sparsesGen['GenSecPeakPrompt']:
                print(f'ERROR: sparse {inputCfg["sparsenameGenPromptSecPeak"]} not found!')
                return None, None
            if not no_List_isMC:
                sparsesGen['GenSecPeakFD'] = inlistData.FindObject(inputCfg['sparsenameGenFDSecPeak'])
            else:
                sparsesGen['GenSecPeakFD'] = indirData.Get(inputCfg['sparsenameGenFDSecPeak'])
                sparsesGen['GenSecPeakFD'].SetName("SparseGenSecPeakFD")
            if not sparsesGen['GenSecPeakFD']:
                print(f'ERROR: sparse {inputCfg["sparsenameGenFDSecPeak"]} not found!')
                return None, None
    infileData.Close()

    return sparses, sparsesGen


def LoadSingleSparseFromTask(infilename, inputCfg, sparsetype='sparsenameBkg'):
    '''
    Method to retrieve single sparse from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file
    - sparse that should be returned

    Returns
    ----------
    - selected sparse from file
    '''
    print('Loading THnSparse from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None
    inlistData = indirData.Get(inputCfg['listname'])
    if not inlistData:
        print(f'List {inputCfg["listname"]} not found!')
        return None
    sparse = inlistData.FindObject(inputCfg[sparsetype])
    if not sparse:
        print(f'ERROR: sparse {inputCfg[sparsetype]} not found!')
        return None

    return sparse


def LoadNormObjFromTask(infilename, inputCfg):
    '''
    Method to retrieve normalisation objects from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - histo with event info and normalisation counter
    '''
    print('Loading norm objects from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None, None
    inlistData = indirData.Get(inputCfg['listname'])
    if not inlistData:
        print(f'List {inputCfg["listname"]} not found!')
        return None, None
    normCounter = indirData.Get(inputCfg['normname'])
    if not normCounter:
        print(f'Norm counter {inputCfg["normname"]} not found!')
        return None, None
    hEv = inlistData.FindObject(inputCfg['histoevname'])
    if not hEv:
        print(f'Histogram {inputCfg["histoevname"]} not found!')
        return None, None

    return hEv, normCounter


def LoadListFromTask(infilename, inputCfg):
    '''
    Method to retrieve list of objects from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - TList of objects from input file
    '''
    print('Loading TList from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None
    inlistData = indirData.Get(inputCfg['listname'])
    if not inlistData:
        print(f'List {inputCfg["listname"]} not found!')
        return None

    return inlistData


def LoadCutObjFromTask(infilename, inputCfg):
    '''
    Method to retrieve D-meson cut object from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - D-meson cut object
    '''
    print('Loading cut object from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None, None
    cutobjname = inputCfg['listname'].replace('coutputDs', 'coutputDsCuts')
    cutobjname = cutobjname.replace('coutputDplus', 'coutputDplusCuts')
    cutobj = indirData.Get(cutobjname)
    if not cutobj:
        print(f'Cut object {cutobjname} not found!')
        return None, None

    return cutobj, cutobjname


def LoadPIDSparses(infilename, inputCfg):
    '''
    Method to retrieve PID sparses from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - sparse of NsigmaTPC and NsigmaTOF variables
    - sparse of NsigmaComb variables
    - dictionary with variable : sparse-axis number
    '''
    print('Loading PID THnSparses from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None, None, None
    inlistData = indirData.Get(inputCfg['listname'])
    if not inlistData:
        print(f'List {inputCfg["listname"]} not found!')
        return None, None, None
    sparsePIDNsigma = inlistData.FindObject('fnSparsePID')
    if not sparsePIDNsigma:
        print('ERROR: sparse fnSparsePID not found!')
        return None, None, None
    sparsePIDNsigmaComb = inlistData.FindObject('fnSparsePIDcomb')
    if not sparsePIDNsigma:
        print('ERROR: sparse fnSparsePIDcomb not found!')
        return None, None, None

    # dictionary of sparse axes with detectors, mass hypothesis, daughter number
    axes = {'TPC': {'Pi': {'0': 2, '1': 6, '2': 10}, 'K': {'0': 3, '1': 7, '2': 11}},
            'TOF': {'Pi': {'0': 4, '1': 8, '2': 12}, 'K': {'0': 5, '1': 9, '2': 13}},
            'Comb': {'Pi': {'0': 2, '1': 4, '2': 6}, 'K': {'0': 3, '1': 5, '2': 7}}}

    return sparsePIDNsigma, sparsePIDNsigmaComb, axes


def LoadSparseFromTaskV2(inputCfg):
    '''
    Method to retrieve sparses from D-meson vn output task file

    Inputs
    ----------
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - list of sparses with reconstructed quantities for all candidates
    '''
    infilename = inputCfg['filename']
    print(f'Loading THnSparses from file {infilename}')
    sparses = []
    infileData = TFile(infilename)

    for dirname, listname in zip(inputCfg['dirname'], inputCfg['listname']):
        indirData = infileData.Get(dirname)
        if not indirData:
            print(f'Directory {dirname} not found!')
            return []
        inlistData = indirData.Get(listname)
        if not inlistData:
            print(f'List {listname} not found!')
            return []
        sparse = inlistData.FindObject(inputCfg['sparsename'])
        if not sparse:
            print('Sparse not found!')
            return []
        sparses.append(sparse)

    return sparses


def LoadListFromTaskV2(infilename, dirname, listname):
    '''
    Method to retrieve list of objects from D-meson vn output task file

    Inputs
    ----------
    - input root file name
    - name of directory in input root file
    - name of list in directory in input root file

    Returns
    ----------
    - TList of objects from input file
    '''
    print('Loading TList from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(dirname)
    if not indirData:
        print(f'Directory {dirname} not found!')
        return None
    inlistData = indirData.Get(listname)
    if not inlistData:
        print(f'List {listname} not found!')
        return None

    return inlistData
