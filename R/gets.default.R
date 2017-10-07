gets.default <-
function(x, ...){
  warning(paste("The 'gets' method does not know how to handle objects of class '",
    class(x), "'", sep=""))
}
