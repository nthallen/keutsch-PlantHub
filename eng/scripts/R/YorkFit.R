#"YorkFit", written by Rick Wehr, 2011, translated into R by Rachel Chang

#Universal routine for finding the best straight line fit
#to data with variable, correlated errors,
#including error and goodness of fit estimates, following Eq. (13) of
#York 2004, American Journal of Physics, which was based in turn on
#York 1969, Earth and Planetary Sciences Letters

YorkFit <- function(X,Y, Xstd, Ystd, Ri=0, b0=0, printCoefs=0, makeLine=0,eps=1e-7){

#X, Y, Xstd, Ystd: waves containing X points, Y points, and their standard deviations
#  WARNING: Xstd and Ystd cannot be zero as this will cause Xw or Yw to be NaN
#  Use a very small value instead.
#Ri: wave containing correlation coefficients for X and Y errors
#b0: rough initial guess for the slope (can be gotten from a standard least-squares fit without errors)
#printCoefs: set equal to 1 to display results in the command window
#makeLine: set equal to 1 to generate a Y wave for the fit line
#Returns a matrix with the intercept and slope plus their uncertainties
  
#If no initial guess for b0 is provided, then just use OLS
  if (b0 == 0) {b0 = lm(Y~X)$coefficients[2]} 

  tol = abs(b0)*eps #the fit will stop iterating when the slope converges to within this value
  
#a,b: final intercept and slope
#a.err, b.err: estimated uncertainties in intercept and slope

  
  # WAVE DEFINITIONS #
  
  Xw = 1/(Xstd^2) #X weights
  Yw = 1/(Ystd^2) #Y weights
  
     
  # ITERATIVE CALCULATION OF SLOPE AND INTERCEPT #
  
  b = b0
  b.diff = tol + 1
  while(b.diff>tol)
{
    b.old = b
    alpha.i = sqrt(Xw*Yw)
    Wi = (Xw*Yw)/((b^2)*Yw + Xw - 2*b*Ri*alpha.i)
    WiX = Wi*X
    WiY = Wi*Y
    sumWiX = sum(WiX, na.rm = TRUE)
    sumWiY = sum(WiY, na.rm = TRUE)
    sumWi = sum(Wi, na.rm = TRUE)
    Xbar = sumWiX/sumWi
    Ybar = sumWiY/sumWi
    Ui = X - Xbar
    Vi = Y - Ybar
    
    Bi = Wi*((Ui/Yw) + (b*Vi/Xw) - (b*Ui+Vi)*Ri/alpha.i)
    wTOPint = Bi*Wi*Vi
    wBOTint = Bi*Wi*Ui
    sumTOP = sum(wTOPint, na.rm=TRUE)
    sumBOT = sum(wBOTint, na.rm=TRUE)
    b = sumTOP/sumBOT
    
    b.diff = abs(b-b.old)
}   
  
  a = Ybar - b*Xbar
  wYorkFitCoefs = c(a,b)
  
  # ERROR CALCULATION #
  
  Xadj = Xbar + Bi
  WiXadj = Wi*Xadj
  sumWiXadj = sum(WiXadj, na.rm=TRUE)
  Xadjbar = sumWiXadj/sumWi
  Uadj = Xadj - Xadjbar
  wErrorTerm = Wi*Uadj*Uadj
  errorSum = sum(wErrorTerm, na.rm=TRUE)
  b.err = sqrt(1/errorSum)
  a.err = sqrt((1/sumWi) + (Xadjbar^2)*(b.err^2))
  wYorkFitErrors = c(a.err,b.err)
  
  # GOODNESS OF FIT CALCULATION #
  lgth = length(X)
  wSint = Wi*(Y - b*X - a)^2
  sumSint = sum(wSint, na.rm=TRUE)
  wYorkGOF = c(sumSint/(lgth-2),sqrt(2/(lgth-2))) #GOF (should equal 1 if assumptions are valid), #standard error in GOF
  
  # OPTIONAL OUTPUTS #
  
  if(printCoefs==1)
{
    print(paste("intercept = ", a, " +/- ", a.err, sep=""))
    print(paste("slope = ", b, " +/- ", b.err, sep=""))
}
  if(makeLine==1)
{
    wYorkFitLine = a + b*X
}
        res=rbind(c(a,a.err),c(b, b.err)) ; dimnames(res)=list(c("a","b"),c("Value","Std.Err"))
  return(res)
}

