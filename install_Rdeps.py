from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
import rpy2.robjects.packages as rpackages

utils = rpackages.importr('utils')
utils.chooseBioCmirror(ind=1)
utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
packnames = ('BiocManager', 'BiocParallel')
names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]

if len(names_to_install) > 0:
    utils.install_packages(StrVector(names_to_install))

BiocManager = importr('BiocManager')
BiocManager.install('SingleR')
BiocManager.install('celldex')
