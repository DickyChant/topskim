from HiggsAnalysis.CombinedLimit.PhysicsModel import *
from HiggsAnalysis.CombinedLimit.NuisanceModifier import *

class TopBtagInMediumModel(PhysicsModel):

    """
    implements a parametric signal scaling based on the b-finding efficiency
    assuming that one b-jet crosses the medium and the other doesn't
    """

    def __init__(self): 
        
        #analysis-specific parameters
        self.A=208
        self.sigma_pp=69e-6
        self.delta_acc=0.00037425
        self.eb_exp=[0.4214,0.415]
        self.eb_exp_uncs={'eb_beff':[0.002,0.002],
                          'eb_udsgeff':[0.0003,0.0003],
                          'eb_jecjer':[0.003,0.003],
                          'eb_quench':[0.035,0.044]}

    def setModelBuilder(self, modelBuilder):

        PhysicsModel.setModelBuilder(self, modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):

        "return 1 for backgrounds, r_ib for signal with i=0,1,2 depending on the bin"

        if not self.DC.isSignal[process]:
            return 1.

        rScale='r_'
        for tag in ['0b','1b','2b']:
            if not tag in bin:
                continue
            rScale+=tag
            break        
        return rScale
                
    
    def setPhysicsOptions(self,physOptions):

        """ parse physics options passed through command line """

        for po in physOptions:
            if 'delta_acc' in po:
                self.delta_acc=float(po.split('=')[1])
                print('Acceptance uncertainty',self.delta_acc)
            if 'sigma_pp' in po:
                self.sigma_pp=float(po.split('=')[1])
                print('Reference pp xsec=',self.sigma_pp)
            if 'A' in po:
                self.A=float(po.split('=')[1])
                print('A=',self.A)
                                                                                                            
    def doParametersOfInterest(self):

        """Create POI and other parameters, and define the POI set."""

        signals=','.join(self.DC.signals)

        #signal strength
        sigmaExp=(self.A**2)*self.sigma_pp
        self.modelBuilder.doVar("sigma[%f,0,%f]"%(sigmaExp,10*sigmaExp))
        self.modelBuilder.doVar("sigma_pp[%f]"%(self.sigma_pp))
        self.modelBuilder.doVar("A[%f]"%(self.A))
        self.modelBuilder.factory_('expr::mu("@0/(@1*@1*@2)",sigma,A,sigma_pp)')
        doAddNuisance(self.DC, (signals, '*', 'accUnc',   'lnN', 1+self.delta_acc) )

        #nuisances for b-tagging
        ebNuisList=[]
        for key in self.eb_exp_uncs:
            nuisName=key+'_nuis'            
            self.modelBuilder.doVar("%s[0,-5,5]"%nuisName)
            ebNuisList.append(nuisName)

        #loop over the jets
        for ij in range(2):
            
            #central value
            self.modelBuilder.doVar("eb%d_cen[%f]"%(ij,self.eb_exp[ij]))            

            #start uncertainties for 
            for key in self.eb_exp_uncs:
                self.modelBuilder.doVar("%s_sigma_%d[%f]"%(key,ij,self.eb_exp_uncs[key][ij]))
            
            #build the final expression to compute the b-finding efficiency
            expr='@0'
            iparam=1
            for key in self.eb_exp_uncs:
                expr += '*(1+%f*@%d)'%(self.eb_exp_uncs[key][ij],iparam)
                iparam+=1
            exprVars='eb%d_cen,'%ij + ','.join(ebNuisList)
            self.modelBuilder.factory_('expr::eb%d("%s",%s)'%(ij,expr,exprVars))
                
        #final scaling expressions
        self.modelBuilder.factory_('expr::r_0b("@0*(1-@1)*(1-@2)",mu,eb0,eb1)')
        self.modelBuilder.factory_('expr::r_1b("@0*(@1+@2-2*@1*@2)",mu,eb0,eb1)')
        self.modelBuilder.factory_('expr::r_2b("@0*@1*@2",mu,eb0,eb1)')

        self.modelBuilder.doSet('POI','sigma')




topBtags = TopBtagInMediumModel()

