# Fixed 19-point grid: lambda_i = 0.25 + 0.75 i/20, i=1,...,19.
# Calibrated with 5 million Gaussian bridges and the independent endpoint
# integrated analytically. Values are stored at full precision.
cv_lambda20_fixed <- function() data.frame(
  m=20L,grid="0.25 + 0.75 i/20, i=1,...,19",
  q0.90_quad=6.255264841259059,q0.95_quad=8.694475491210449,q0.99_quad=14.414683004271287,
  q0.90_range=2.3846824930497283,q0.95_range=3.1993661030719607,q0.99_range=5.023448797423727)
get_cv_lambda20 <- function(alpha=.05,which=c("quad","range")) {
  which<-match.arg(which)
  if(!alpha %in% c(.1,.05,.01))stop("Unsupported alpha")
  cv_lambda20_fixed()[[paste0("q",sprintf("%.2f",1-alpha),"_",which)]]
}
