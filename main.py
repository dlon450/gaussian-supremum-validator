from formulations import *

class Validators:
    def __init__(self, formulation_name) -> None:
        self.formulation = {'dro': CCP_DRO_moment, 'ro_ellipsoid': CCP_RO_ellipsoid, 'saa': CCP_SAA, 'so': CCP_SO}[formulation_name]
    def cross_validation(self):
        pass
    def bootstrapping(self):
        pass
    def sectioning(self):
        pass
    def gaussian(self):
        pass