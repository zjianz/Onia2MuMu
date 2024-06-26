from FWCore.GuiBrowsers.ConfigToolBase import *

from PhysicsTools.PatAlgos.tools.helpers import *

class RemoveAllPATObjectsBut(ConfigToolBase):

    """ Remove all PAT objects from the default sequence but a specific one
    """
    _label='removeAllPATObjectsBut'
    _defaultParameters=dicttypes.SortedKeysDict()
    def __init__(self):
        ConfigToolBase.__init__(self)
        self.addParameter(self._defaultParameters,'names',self._defaultValue, "list of collection names; supported are 'Photons', 'Electrons', 'Muons', 'Taus', 'Jets', 'METs'", Type=list, allowedValues=['Photons', 'Electrons', 'Muons', 'Taus', 'Jets', 'METs'])
        self.addParameter(self._defaultParameters,'outputModules',['out'], "names of all output modules specified to be adapted (default is ['out'])")
        self._parameters=copy.deepcopy(self._defaultParameters)
        self._comment = ""

    def getDefaultParameters(self):
        return self._defaultParameters

    def __call__(self,process,
                 names           = None,
                 outputInProcess = None,
                 outputModules   = None) :
        ## stop processing if 'outputInProcess' exists and show the new alternative
        if  not outputInProcess is None:
            deprecatedOptionOutputInProcess(self)
        if  names is None:
            names=self._defaultParameters['names'].value
        if  outputModules is None:
            outputModules=self._defaultParameters['outputModules'].value
        self.setParameter('names',names)
        self.setParameter('outputModules',outputModules)
        self.apply(process)

    def toolCode(self, process):
        names=self._parameters['names'].value
        outputModules=self._parameters['outputModules'].value

        removeTheseObjectCollections = ['Photons', 'Electrons', 'Muons', 'Taus', 'Jets', 'METs']
        for obj in range(len(names)):
            removeTheseObjectCollections.remove(names[obj])
        removeSpecificPATObjects(process, removeTheseObjectCollections, outputModules = outputModules)

removeAllPATObjectsBut=RemoveAllPATObjectsBut()


class RemoveSpecificPATObjects(ConfigToolBase):

    """ Remove a specific PAT object from the default sequence
    """
    _label='removeSpecificPATObjects'
    _defaultParameters=dicttypes.SortedKeysDict()
    def __init__(self):
        ConfigToolBase.__init__(self)
        self.addParameter(self._defaultParameters,'names',self._defaultValue, "list of collection names; supported are 'Photons', 'Electrons', 'Muons', 'Taus', 'Jets', 'METs'", Type=list, allowedValues=['Photons', 'Electrons', 'Muons', 'Taus', 'Jets', 'METs'])
        self.addParameter(self._defaultParameters,'outputModules',['out'], "names of all output modules specified to be adapted (default is ['out'])")
        self.addParameter(self._defaultParameters,'postfix',"", "postfix of default sequence")
        self._parameters=copy.deepcopy(self._defaultParameters)
        self._comment = ""

    def getDefaultParameters(self):
        return self._defaultParameters

    def __call__(self,process,
                 names           = None,
                 outputInProcess = None,
                 postfix         = None,
                 outputModules   = None) :
        ## stop processing if 'outputInProcess' exists and show the new alternative
        if  not outputInProcess is None:
            deprecatedOptionOutputInProcess(self)
        if  names is None:
            names=self._defaultParameters['names'].value
        if  outputModules is None:
            outputModules=self._defaultParameters['outputModules'].value
        if postfix  is None:
            postfix=self._defaultParameters['postfix'].value
        self.setParameter('names',names)
        self.setParameter('outputModules',outputModules)
        self.setParameter('postfix',postfix)
        self.apply(process)

    def toolCode(self, process):
        names=self._parameters['names'].value
        outputModules=self._parameters['outputModules'].value
        postfix=self._parameters['postfix'].value

        ## remove pre object production steps from the default sequence
        for obj in range(len(names)):
            if( names[obj] == 'Photons' ):
                removeIfInSequence(process, 'patPhotonIsolation', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'photonMatch', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patPhotons', "patDefaultSequence", postfix)
            if( names[obj] == 'Electrons' ):
                removeIfInSequence(process, 'patElectronId', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patElectronIsolation', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'electronMatch', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patElectrons', "patDefaultSequence", postfix)
            if( names[obj] == 'Muons' ):
                removeIfInSequence(process, 'muonMatch', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patMuons', "patDefaultSequence", postfix)
            if( names[obj] == 'Taus' ):
                removeIfInSequence(process, 'patPFCandidateIsoDepositSelection', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patPFTauIsolation', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'tauMatch', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'tauGenJets', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'tauGenJetsSelectorAllHadrons', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'tauGenJetMatch', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patTaus', "patDefaultSequence", postfix)
            if( names[obj] == 'Jets' ):
                removeIfInSequence(process, 'patJetCharge', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patJetCorrections', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patJetPartonMatch', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patJetGenJetMatch', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patJetFlavourIdLegacy', "patDefaultSequence", postfix)
                removeIfInSequence(process, 'patJetFlavourId', "patDefaultSequence", postfix)
            if( names[obj] == 'METs' ):
                removeIfInSequence(process, 'patMETCorrections', "patDefaultSequence", postfix)

            ## remove object production steps from the default sequence
            if( names[obj] == 'METs' ):
                process.patDefaultSequence.remove( getattr(process, 'pat'+names[obj]) )
            else:
                if( names[obj] == 'Jets' ):
                    applyPostfix(process,"patDefaultSequence",postfix).remove(
                        getattr(process, jetCollectionString()+postfix) )
                    applyPostfix(process,"patDefaultSequence",postfix).remove(
                        getattr(process, jetCollectionString('selected')+postfix) )
                    applyPostfix(process,"patDefaultSequence",postfix).remove(
                        getattr(process, jetCollectionString('count')+postfix) )
                else:
                    applyPostfix(process,"patDefaultSequence",postfix).remove(
                        getattr(process, 'pat'+names[obj]+postfix) )
                    applyPostfix(process,"patDefaultSequence",postfix).remove(
                        getattr(process, 'selectedPat'+names[obj]+postfix) )
                    applyPostfix(process,"patDefaultSequence",postfix).remove(
                        getattr(process, 'countPat'+names[obj]+postfix) )
            ## in the case of leptons, the lepton counter must be modified as well
            if( names[obj] == 'Electrons' ):
                print('removed from lepton counter: electrons')
                applyPostfix(process,"countPatLeptons",postfix).countElectrons = False
            elif( names[obj] == 'Muons' ):
                print('removed from lepton counter: muons')
                applyPostfix(process,"countPatLeptons",postfix).countMuons = False
            elif( names[obj] == 'Taus' ):
                print('removed from lepton counter: taus')
                applyPostfix(process,"countPatLeptons",postfix).countTaus = False
            ## remove from summary
            if( names[obj] == 'METs' ):
                applyPostfix(process,"patCandidateSummary",postfix).candidates.remove(
                    cms.InputTag('pat'+names[obj]+postfix) )
            else:
                if( names[obj] == 'Jets' ):
                    applyPostfix(process,"patCandidateSummary",postfix).candidates.remove(
                        cms.InputTag(jetCollectionString()+postfix) )
                    applyPostfix(process,"selectedPatCandidateSummary",postfix).candidates.remove(
                        cms.InputTag(jetCollectionString('selected')+postfix) )
                    applyPostfix(process,"cleanPatCandidateSummary",postfix).candidates.remove(
                        cms.InputTag(jetCollectionString('clean')+postfix) )
                else:
                    ## check whether module is in sequence or not
                    result = [ m.label()[:-len(postfix)] for m in listModules( getattr(process,"patDefaultSequence"+postfix))]
                    result.extend([ m.label()[:-len(postfix)] for m in listSequences( getattr(process,"patDefaultSequence"+postfix))]  )
                    if applyPostfix(process,"patCandidateSummary",postfix) in result :
                        applyPostfix(process,"patCandidateSummary",postfix).candidates.remove(
                            cms.InputTag('pat'+names[obj]+postfix) )
                    if applyPostfix(process,"selectedPatCandidateSummary",postfix) in result :
                        applyPostfix(process,"selectedPatCandidateSummary",postfix).candidates.remove(
                            cms.InputTag('selectedPat'+names[obj]+postfix) )
                    if applyPostfix(process,"cleanPatCandidateSummary",postfix) in result :
                        applyPostfix(process,"cleanPatCandidateSummary",postfix).candidates.remove(
                            cms.InputTag('cleanPat'+names[obj]+postfix) )
        ## remove cleaning for the moment; in principle only the removed object
        ## could be taken out of the checkOverlaps PSet
        if len(outputModules) > 0:
            print("---------------------------------------------------------------------")
            print("INFO   : some objects have been removed from the sequence. Switching ")
            print("         off PAT cross collection cleaning, as it might be of limited")
            print("         sense now. If you still want to keep object collection cross")
            print("         cleaning within PAT you need to run and configure it by hand")
            removeCleaning(process,outputModules=outputModules,postfix=postfix)

removeSpecificPATObjects=RemoveSpecificPATObjects()

class RemoveCleaning(ConfigToolBase):

    """ remove PAT cleaning from the default sequence:
    """
    _label='removeCleaning'
    _defaultParameters=dicttypes.SortedKeysDict()
    def __init__(self):
        ConfigToolBase.__init__(self)
        self.addParameter(self._defaultParameters,'outputModules',['out'], "names of all output modules specified to be adapted (default is ['out'])")
        self.addParameter(self._defaultParameters,'postfix',"", "postfix of default sequence")
        self._parameters=copy.deepcopy(self._defaultParameters)
        self._comment = ""

    def getDefaultParameters(self):
        return self._defaultParameters

    def __call__(self,process,
                 outputInProcess = None,
                 postfix         = None,
                 outputModules   = None) :
        ## stop processing if 'outputInProcess' exists and show the new alternative
        if  not outputInProcess is None:
            deprecatedOptionOutputInProcess(self)
        if  outputModules is None:
            outputModules=self._defaultParameters['outputModules'].value
        if postfix  is None:
            postfix=self._defaultParameters['postfix'].value

        self.setParameter('outputModules',outputModules)
        self.setParameter('postfix',postfix)

        self.apply(process)

    def toolCode(self, process):
        outputModules=self._parameters['outputModules'].value
        postfix=self._parameters['postfix'].value

        ## adapt single object counters
        for m in listModules(applyPostfix(process,"countPatCandidates",postfix)):
            if hasattr(m, 'src'): m.src = m.src.value().replace('cleanPat','selectedPat')

        ## adapt lepton counter
        countLept = applyPostfix(process,"countPatLeptons",postfix)
        countLept.electronSource = countLept.electronSource.value().replace('cleanPat','selectedPat')
        countLept.muonSource = countLept.muonSource.value().replace('cleanPat','selectedPat')
        countLept.tauSource = countLept.tauSource.value().replace('cleanPat','selectedPat')
        for m in getattr(process, "cleanPatCandidates").moduleNames():
            getattr(process, "patDefaultSequence"+postfix).remove(
                applyPostfix(process,m,postfix)
                )
        if len(outputModules) > 0:
            print("------------------------------------------------------------")
            print("INFO   : cleaning has been removed. Switching output from")
            print("         clean PAT candidates to selected PAT candidates.")
            ## add selected pat objects to the pat output
            from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
            for outMod in outputModules:
                if hasattr(process,outMod):
                    getattr(process,outMod).outputCommands = patEventContentNoCleaning
                else:
                    raise KeyError("process has no OutModule named", outMod)

removeCleaning=RemoveCleaning()
