## Script to create counterbalanced task order using latin square for COLA RR
## Created Z. Woodhead 24/11/2020

ntasks <- 6
nsubj <- 48 # Must be a multiple of ntasks
seeds <- c(100, 200, 300, 400, 500) # One seed for each of 5 site
whichsite <- as.integer(readline(prompt="Which site? 1=Bangor, 2=Lancaster, 3=Lincoln, 4=Oxford, 5=UCL: "))           
myseed <- seeds[whichsite]

## Latinsquare function from http://www.cookbook-r.com/Tools_for_experiments/Generating_counterbalanced_orders/
## - len is the size of the latin square
## - reps is the number of repetitions - how many Latin squares to generate
## - seed is a random seed that can be used to generate repeatable sequences
## - returnstrings tells it to return a vector of char strings for each square,
##    instead of a big matrix. This option is only really used for checking the
##    randomness of the squares.
latinsquare <- function(len, reps=1, seed=NA, returnstrings=FALSE) {
  
  # Save the old random seed and use the new one, if present
  if (!is.na(seed)) {
    if (exists(".Random.seed"))  { saved.seed <- .Random.seed }
    else                         { saved.seed <- NA }
    set.seed(seed)
  }
  
  # This matrix will contain all the individual squares
  allsq <- matrix(nrow=reps*len, ncol=len)
  
  # Store a string id of each square if requested
  if (returnstrings) {  squareid <- vector(mode = "character", length = reps) }
  
  # Get a random element from a vector (the built-in sample function annoyingly
  #   has different behavior if there's only one element in x)
  sample1 <- function(x) {
    if (length(x)==1) { return(x) }
    else              { return(sample(x,1)) }
  }
  
  # Generate each of n individual squares
  for (n in 1:reps) {
    
    # Generate an empty square
    sq <- matrix(nrow=len, ncol=len) 
    
    # If we fill the square sequentially from top left, some latin squares
    # are more probable than others.  So we have to do it random order,
    # all over the square.
    # The rough procedure is:
    # - randomly select a cell that is currently NA (call it the target cell)
    # - find all the NA cells sharing the same row or column as the target
    # - fill the target cell
    # - fill the other cells sharing the row/col
    # - If it ever is impossible to fill a cell because all the numbers
    #    are already used, then quit and start over with a new square.
    # In short, it picks a random empty cell, fills it, then fills in the 
    # other empty cells in the "cross" in random order. If we went totally randomly
    # (without the cross), the failure rate is much higher.
    while (any(is.na(sq))) {
      
      # Pick a random cell which is currently NA
      k <- sample1(which(is.na(sq)))
      
      i <- (k-1) %% len +1       # Get the row num
      j <- floor((k-1) / len) +1 # Get the col num
      
      # Find the other NA cells in the "cross" centered at i,j
      sqrow <- sq[i,]
      sqcol <- sq[,j]
      
      # A matrix of coordinates of all the NA cells in the cross
      openCell <-rbind( cbind(which(is.na(sqcol)), j),
                        cbind(i, which(is.na(sqrow))))
      # Randomize fill order
      openCell <- openCell[sample(nrow(openCell)),]
      
      # Put center cell at top of list, so that it gets filled first
      openCell <- rbind(c(i,j), openCell)
      # There will now be three entries for the center cell, so remove duplicated entries
      # Need to make sure it's a matrix -- otherwise, if there's just 
      # one row, it turns into a vector, which causes problems
      openCell <- matrix(openCell[!duplicated(openCell),], ncol=2)
      
      # Fill in the center of the cross, then the other open spaces in the cross
      for (c in 1:nrow(openCell)) {
        # The current cell to fill
        ci <- openCell[c,1]
        cj <- openCell[c,2]
        # Get the numbers that are unused in the "cross" centered on i,j
        freeNum <- which(!(1:len %in% c(sq[ci,], sq[,cj])))
        
        # Fill in this location on the square
        if (length(freeNum)>0) { sq[ci,cj] <- sample1(freeNum) }
        else  {
          # Failed attempt - no available numbers
          # Re-generate empty square
          sq <- matrix(nrow=len, ncol=len)
          
          # Break out of loop
          break;
        }
      }
    }
    
    # Store the individual square into the matrix containing all squares
    allsqrows <- ((n-1)*len) + 1:len
    allsq[allsqrows,] <- sq
    
    # Store a string representation of the square if requested. Each unique
    # square has a unique string.
    if (returnstrings) { squareid[n] <- paste(sq, collapse="") }
    
  }
  
  # Restore the old random seed, if present
  if (!is.na(seed) && !is.na(saved.seed)) { .Random.seed <- saved.seed }
  
  if (returnstrings) { return(squareid) }
  else               { return(allsq) }
}


## Create task orders 
taskorder <- latinsquare(ntasks, reps=(nsubj/ntasks), seed = myseed)
taskorder_letters <- matrix(data=NA, nrow=(dim(taskorder)[1]), ncol=(dim(taskorder)[2]))
for (t in 1:ntasks){
  taskorder_letters[, t] <- LETTERS[taskorder[ , t]]
}

taskorder_letters <- as.data.frame((taskorder_letters))
colnames(taskorder_letters) <- c('T1','T2','T3','T4','T5','T6')
taskorder_letters$ID <- seq(from=(myseed+1), to=(myseed+nsubj))
taskorder_letters$list <- with(taskorder_letters, paste0(T1,T2,T3,T4,T5,T6))
sites <- c('Bangor','Lancaster','Lincoln','Oxford','UCL')

write.csv(taskorder_letters, paste0('taskorder_', sites[whichsite], '.csv'))
