# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.
import _B03App
def _swig_setattr(self,class_type,name,value):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    self.__dict__[name] = value

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class G4RunManager(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4RunManager, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4RunManager, name)
    def BeamOn(*args): return apply(_B03App.G4RunManager_BeamOn,args)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C G4RunManager instance at %s>" % (self.this,)

class G4RunManagerPtr(G4RunManager):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4RunManager
_B03App.G4RunManager_swigregister(G4RunManagerPtr)

class G4VIStore(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4VIStore, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4VIStore, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C G4VIStore instance at %s>" % (self.this,)

class G4VIStorePtr(G4VIStore):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4VIStore
_B03App.G4VIStore_swigregister(G4VIStorePtr)

class G4VPhysicalVolume(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4VPhysicalVolume, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4VPhysicalVolume, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C G4VPhysicalVolume instance at %s>" % (self.this,)

class G4VPhysicalVolumePtr(G4VPhysicalVolume):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4VPhysicalVolume
_B03App.G4VPhysicalVolume_swigregister(G4VPhysicalVolumePtr)

class G4GeometryCell(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4GeometryCell, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4GeometryCell, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_G4GeometryCell,args)
        self.thisown = 1
    def __del__(self, destroy= _B03App.delete_G4GeometryCell):
        try:
            if self.thisown: destroy(self)
        except: pass
    def GetPhysicalVolume(*args): return apply(_B03App.G4GeometryCell_GetPhysicalVolume,args)
    def GetReplicaNumber(*args): return apply(_B03App.G4GeometryCell_GetReplicaNumber,args)
    def __repr__(self):
        return "<C G4GeometryCell instance at %s>" % (self.this,)

class G4GeometryCellPtr(G4GeometryCell):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4GeometryCell
_B03App.G4GeometryCell_swigregister(G4GeometryCellPtr)

class G4IStore(G4VIStore):
    __swig_setmethods__ = {}
    for _s in [G4VIStore]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4IStore, name, value)
    __swig_getmethods__ = {}
    for _s in [G4VIStore]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, G4IStore, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_G4IStore,args)
        self.thisown = 1
    def __del__(self, destroy= _B03App.delete_G4IStore):
        try:
            if self.thisown: destroy(self)
        except: pass
    def AddImportanceGeometryCell(*args): return apply(_B03App.G4IStore_AddImportanceGeometryCell,args)
    def ChangeImportance(*args): return apply(_B03App.G4IStore_ChangeImportance,args)
    def GetImportance(*args): return apply(_B03App.G4IStore_GetImportance,args)
    def IsKnown(*args): return apply(_B03App.G4IStore_IsKnown,args)
    def GetWorldVolume(*args): return apply(_B03App.G4IStore_GetWorldVolume,args)
    def __repr__(self):
        return "<C G4IStore instance at %s>" % (self.this,)

class G4IStorePtr(G4IStore):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4IStore
_B03App.G4IStore_swigregister(G4IStorePtr)

class G4CellScoreValues(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4CellScoreValues, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4CellScoreValues, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_G4CellScoreValues,args)
        self.thisown = 1
    __swig_setmethods__["fSumSL"] = _B03App.G4CellScoreValues_fSumSL_set
    __swig_getmethods__["fSumSL"] = _B03App.G4CellScoreValues_fSumSL_get
    if _newclass:fSumSL = property(_B03App.G4CellScoreValues_fSumSL_get,_B03App.G4CellScoreValues_fSumSL_set)
    __swig_setmethods__["fSumSLW"] = _B03App.G4CellScoreValues_fSumSLW_set
    __swig_getmethods__["fSumSLW"] = _B03App.G4CellScoreValues_fSumSLW_get
    if _newclass:fSumSLW = property(_B03App.G4CellScoreValues_fSumSLW_get,_B03App.G4CellScoreValues_fSumSLW_set)
    __swig_setmethods__["fSumSLW_v"] = _B03App.G4CellScoreValues_fSumSLW_v_set
    __swig_getmethods__["fSumSLW_v"] = _B03App.G4CellScoreValues_fSumSLW_v_get
    if _newclass:fSumSLW_v = property(_B03App.G4CellScoreValues_fSumSLW_v_get,_B03App.G4CellScoreValues_fSumSLW_v_set)
    __swig_setmethods__["fSumSLWE"] = _B03App.G4CellScoreValues_fSumSLWE_set
    __swig_getmethods__["fSumSLWE"] = _B03App.G4CellScoreValues_fSumSLWE_get
    if _newclass:fSumSLWE = property(_B03App.G4CellScoreValues_fSumSLWE_get,_B03App.G4CellScoreValues_fSumSLWE_set)
    __swig_setmethods__["fSumSLWE_v"] = _B03App.G4CellScoreValues_fSumSLWE_v_set
    __swig_getmethods__["fSumSLWE_v"] = _B03App.G4CellScoreValues_fSumSLWE_v_get
    if _newclass:fSumSLWE_v = property(_B03App.G4CellScoreValues_fSumSLWE_v_get,_B03App.G4CellScoreValues_fSumSLWE_v_set)
    __swig_setmethods__["fSumTracksEntering"] = _B03App.G4CellScoreValues_fSumTracksEntering_set
    __swig_getmethods__["fSumTracksEntering"] = _B03App.G4CellScoreValues_fSumTracksEntering_get
    if _newclass:fSumTracksEntering = property(_B03App.G4CellScoreValues_fSumTracksEntering_get,_B03App.G4CellScoreValues_fSumTracksEntering_set)
    __swig_setmethods__["fSumPopulation"] = _B03App.G4CellScoreValues_fSumPopulation_set
    __swig_getmethods__["fSumPopulation"] = _B03App.G4CellScoreValues_fSumPopulation_get
    if _newclass:fSumPopulation = property(_B03App.G4CellScoreValues_fSumPopulation_get,_B03App.G4CellScoreValues_fSumPopulation_set)
    __swig_setmethods__["fSumCollisions"] = _B03App.G4CellScoreValues_fSumCollisions_set
    __swig_getmethods__["fSumCollisions"] = _B03App.G4CellScoreValues_fSumCollisions_get
    if _newclass:fSumCollisions = property(_B03App.G4CellScoreValues_fSumCollisions_get,_B03App.G4CellScoreValues_fSumCollisions_set)
    __swig_setmethods__["fSumCollisionsWeight"] = _B03App.G4CellScoreValues_fSumCollisionsWeight_set
    __swig_getmethods__["fSumCollisionsWeight"] = _B03App.G4CellScoreValues_fSumCollisionsWeight_get
    if _newclass:fSumCollisionsWeight = property(_B03App.G4CellScoreValues_fSumCollisionsWeight_get,_B03App.G4CellScoreValues_fSumCollisionsWeight_set)
    __swig_setmethods__["fNumberWeightedEnergy"] = _B03App.G4CellScoreValues_fNumberWeightedEnergy_set
    __swig_getmethods__["fNumberWeightedEnergy"] = _B03App.G4CellScoreValues_fNumberWeightedEnergy_get
    if _newclass:fNumberWeightedEnergy = property(_B03App.G4CellScoreValues_fNumberWeightedEnergy_get,_B03App.G4CellScoreValues_fNumberWeightedEnergy_set)
    __swig_setmethods__["fFluxWeightedEnergy"] = _B03App.G4CellScoreValues_fFluxWeightedEnergy_set
    __swig_getmethods__["fFluxWeightedEnergy"] = _B03App.G4CellScoreValues_fFluxWeightedEnergy_get
    if _newclass:fFluxWeightedEnergy = property(_B03App.G4CellScoreValues_fFluxWeightedEnergy_get,_B03App.G4CellScoreValues_fFluxWeightedEnergy_set)
    __swig_setmethods__["fAverageTrackWeight"] = _B03App.G4CellScoreValues_fAverageTrackWeight_set
    __swig_getmethods__["fAverageTrackWeight"] = _B03App.G4CellScoreValues_fAverageTrackWeight_get
    if _newclass:fAverageTrackWeight = property(_B03App.G4CellScoreValues_fAverageTrackWeight_get,_B03App.G4CellScoreValues_fAverageTrackWeight_set)
    __swig_setmethods__["fImportance"] = _B03App.G4CellScoreValues_fImportance_set
    __swig_getmethods__["fImportance"] = _B03App.G4CellScoreValues_fImportance_get
    if _newclass:fImportance = property(_B03App.G4CellScoreValues_fImportance_get,_B03App.G4CellScoreValues_fImportance_set)
    def __repr__(self):
        return "<C G4CellScoreValues instance at %s>" % (self.this,)

class G4CellScoreValuesPtr(G4CellScoreValues):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4CellScoreValues
_B03App.G4CellScoreValues_swigregister(G4CellScoreValuesPtr)

class G4CellScoreComposer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4CellScoreComposer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4CellScoreComposer, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_G4CellScoreComposer,args)
        self.thisown = 1
    def __del__(self, destroy= _B03App.delete_G4CellScoreComposer):
        try:
            if self.thisown: destroy(self)
        except: pass
    def GetStandardCellScoreValues(*args): return apply(_B03App.G4CellScoreComposer_GetStandardCellScoreValues,args)
    def __repr__(self):
        return "<C G4CellScoreComposer instance at %s>" % (self.this,)

class G4CellScoreComposerPtr(G4CellScoreComposer):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4CellScoreComposer
_B03App.G4CellScoreComposer_swigregister(G4CellScoreComposerPtr)

class G4CellScorer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4CellScorer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4CellScorer, name)
    def GetCellScoreComposer(*args): return apply(_B03App.G4CellScorer_GetCellScoreComposer,args)
    def GetCellScoreValues(*args): return apply(_B03App.G4CellScorer_GetCellScoreValues,args)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C G4CellScorer instance at %s>" % (self.this,)

class G4CellScorerPtr(G4CellScorer):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4CellScorer
_B03App.G4CellScorer_swigregister(G4CellScorerPtr)

class G4VCellScorerStore(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4VCellScorerStore, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4VCellScorerStore, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C G4VCellScorerStore instance at %s>" % (self.this,)

class G4VCellScorerStorePtr(G4VCellScorerStore):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4VCellScorerStore
_B03App.G4VCellScorerStore_swigregister(G4VCellScorerStorePtr)

class G4CellScorerStore(G4VCellScorerStore):
    __swig_setmethods__ = {}
    for _s in [G4VCellScorerStore]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4CellScorerStore, name, value)
    __swig_getmethods__ = {}
    for _s in [G4VCellScorerStore]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, G4CellScorerStore, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_G4CellScorerStore,args)
        self.thisown = 1
    def __del__(self, destroy= _B03App.delete_G4CellScorerStore):
        try:
            if self.thisown: destroy(self)
        except: pass
    def SetAutoScorerCreate(*args): return apply(_B03App.G4CellScorerStore_SetAutoScorerCreate,args)
    def AddCellScorer(*args): return apply(_B03App.G4CellScorerStore_AddCellScorer,args)
    def GetMapGeometryCellCellScorer(*args): return apply(_B03App.G4CellScorerStore_GetMapGeometryCellCellScorer,args)
    def __repr__(self):
        return "<C G4CellScorerStore instance at %s>" % (self.this,)

class G4CellScorerStorePtr(G4CellScorerStore):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4CellScorerStore
_B03App.G4CellScorerStore_swigregister(G4CellScorerStorePtr)

class G4VScorer(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4VScorer, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4VScorer, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C G4VScorer instance at %s>" % (self.this,)

class G4VScorerPtr(G4VScorer):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4VScorer
_B03App.G4VScorer_swigregister(G4VScorerPtr)

class G4CellStoreScorer(G4VScorer):
    __swig_setmethods__ = {}
    for _s in [G4VScorer]: __swig_setmethods__.update(_s.__swig_setmethods__)
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4CellStoreScorer, name, value)
    __swig_getmethods__ = {}
    for _s in [G4VScorer]: __swig_getmethods__.update(_s.__swig_getmethods__)
    __getattr__ = lambda self, name: _swig_getattr(self, G4CellStoreScorer, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_G4CellStoreScorer,args)
        self.thisown = 1
    def __repr__(self):
        return "<C G4CellStoreScorer instance at %s>" % (self.this,)

class G4CellStoreScorerPtr(G4CellStoreScorer):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4CellStoreScorer
_B03App.G4CellStoreScorer_swigregister(G4CellStoreScorerPtr)

class G4VImportanceAlgorithm(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4VImportanceAlgorithm, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4VImportanceAlgorithm, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C G4VImportanceAlgorithm instance at %s>" % (self.this,)

class G4VImportanceAlgorithmPtr(G4VImportanceAlgorithm):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4VImportanceAlgorithm
_B03App.G4VImportanceAlgorithm_swigregister(G4VImportanceAlgorithmPtr)

class G4ParallelGeometrySampler(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4ParallelGeometrySampler, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4ParallelGeometrySampler, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_G4ParallelGeometrySampler,args)
        self.thisown = 1
    def __del__(self, destroy= _B03App.delete_G4ParallelGeometrySampler):
        try:
            if self.thisown: destroy(self)
        except: pass
    def PrepareScoring(*args): return apply(_B03App.G4ParallelGeometrySampler_PrepareScoring,args)
    def PrepareImportanceSampling(*args): return apply(_B03App.G4ParallelGeometrySampler_PrepareImportanceSampling,args)
    def PrepareWeightRoulett(*args): return apply(_B03App.G4ParallelGeometrySampler_PrepareWeightRoulett,args)
    def Configure(*args): return apply(_B03App.G4ParallelGeometrySampler_Configure,args)
    def ClearSampling(*args): return apply(_B03App.G4ParallelGeometrySampler_ClearSampling,args)
    def IsConfigured(*args): return apply(_B03App.G4ParallelGeometrySampler_IsConfigured,args)
    def __repr__(self):
        return "<C G4ParallelGeometrySampler instance at %s>" % (self.this,)

class G4ParallelGeometrySamplerPtr(G4ParallelGeometrySampler):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4ParallelGeometrySampler
_B03App.G4ParallelGeometrySampler_swigregister(G4ParallelGeometrySamplerPtr)

class G4ScoreTable(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, G4ScoreTable, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, G4ScoreTable, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_G4ScoreTable,args)
        self.thisown = 1
    def __del__(self, destroy= _B03App.delete_G4ScoreTable):
        try:
            if self.thisown: destroy(self)
        except: pass
    def Print(*args): return apply(_B03App.G4ScoreTable_Print,args)
    def Write(*args): return apply(_B03App.G4ScoreTable_Write,args)
    def __repr__(self):
        return "<C G4ScoreTable instance at %s>" % (self.this,)

class G4ScoreTablePtr(G4ScoreTable):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = G4ScoreTable
_B03App.G4ScoreTable_swigregister(G4ScoreTablePtr)

class B03ImportanceDetectorConstruction(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, B03ImportanceDetectorConstruction, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, B03ImportanceDetectorConstruction, name)
    def __init__(self,*args):
        self.this = apply(_B03App.new_B03ImportanceDetectorConstruction,args)
        self.thisown = 1
    def __del__(self, destroy= _B03App.delete_B03ImportanceDetectorConstruction):
        try:
            if self.thisown: destroy(self)
        except: pass
    def GetPhysicalVolumeByName(*args): return apply(_B03App.B03ImportanceDetectorConstruction_GetPhysicalVolumeByName,args)
    def GetWorldVolume(*args): return apply(_B03App.B03ImportanceDetectorConstruction_GetWorldVolume,args)
    def ListPhysNames(*args): return apply(_B03App.B03ImportanceDetectorConstruction_ListPhysNames,args)
    def GetGeometryCell(*args): return apply(_B03App.B03ImportanceDetectorConstruction_GetGeometryCell,args)
    def __repr__(self):
        return "<C B03ImportanceDetectorConstruction instance at %s>" % (self.this,)

class B03ImportanceDetectorConstructionPtr(B03ImportanceDetectorConstruction):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = B03ImportanceDetectorConstruction
_B03App.B03ImportanceDetectorConstruction_swigregister(B03ImportanceDetectorConstructionPtr)

class B03AppBase(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, B03AppBase, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, B03AppBase, name)
    def __del__(self, destroy= _B03App.delete_B03AppBase):
        try:
            if self.thisown: destroy(self)
        except: pass
    __swig_getmethods__["GetB03AppBase"] = lambda x: _B03App.B03AppBase_GetB03AppBase
    if _newclass:GetB03AppBase = staticmethod(_B03App.B03AppBase_GetB03AppBase)
    def GetRunManager(*args): return apply(_B03App.B03AppBase_GetRunManager,args)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<C B03AppBase instance at %s>" % (self.this,)

class B03AppBasePtr(B03AppBase):
    def __init__(self,this):
        self.this = this
        if not hasattr(self,"thisown"): self.thisown = 0
        self.__class__ = B03AppBase
_B03App.B03AppBase_swigregister(B03AppBasePtr)
B03AppBase_GetB03AppBase = _B03App.B03AppBase_GetB03AppBase



