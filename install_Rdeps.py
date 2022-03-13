from rpy2.robjects.packages import importr

devtools = importr('devtools')
base = importr('base')
base.options(timeout=300)  # cisTopic download takes some time
devtools.install_github("aertslab/cisTopic", upgrade=False)
devtools.install_github("CamaraLab/STvEA", force=True, upgrade=False)
