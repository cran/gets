as.lm <-
function(object)
{

  ##what kind of class?:
  objectClass <- class(object)
  classOK <-
    ifelse( objectClass %in% c("arx","gets","isat"), TRUE, FALSE)

  ##class not OK:
  if(!classOK){
    stop("'object' must be of class 'arx', 'gets' or 'isat'")
  }

  ##class OK:
  if(classOK){
    y <- object$aux$y
    x <- object$aux$mX
    colnames(x) <- object$aux$mXnames
    yx <- data.frame(y, x)
    result <- lm(formula = y ~ . - 1, data = yx)
  }
  
  ##return result:
  return(result)
  
}
