## Function: verifyItems
#
# Use this function to quickly generate a T/F vector to annotate or subset a
# given list of items. The function will display each item and ask the user
# to indicate T or F (lowercase is also okay). Invalid user responses will
# trigger a repeat of the question.
#
# list: User provided list of items
# range: User provided subset of list to verify. Default is entire list.
#
# Return: A logical vector the same length as range
#
# Examples:  
#
# my.list <- c("apple","orange","tomato","carrot")
# fruits <- verifyItems(my.list)
# my.list[fruits] # will subset for verified items
#
# my.list <- c("apple","orange","tomato","carrot")
# fruits <- verifyItems(my.list, 1:3) #only review first 3 items
# my.list[1:3][fruits] # will subset for verified items over same range
#
verifyItems <- function(list = NULL, range = NULL){
  if(is.null(list))
    stop("Must provide a list")

  if(is.null(range))
    range <- 1:length(list)
  
  if(!is.vector(list[range]))
    stop("Range must subset provided list.")
  
  sapply(list[range], function(x){
    print(x)
    r<-NULL
    while(!is.logical(r) ){
      r<-toupper(readline("Valid? (T/F): "))
      if (r %in% c("T","F"))
        r <- as.logical(r)
    }
    r
  })
}
