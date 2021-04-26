get_max = function(x) {
  idx = floor(length(x)/2)
  done = FALSE
  i = floor(length(x)/2)
  while (!done) {
    if (x[i] >= x[i + 1] && x[i] >= x[i - 1]){
      return(i)
    } else if (x[i] < x[i + 1]) {
      i = i + 1
    } else {
      i = i - 1
    }
  }
}

x <- seq(from=0,to=1e8,by=1)
start_time <- Sys.time()
y1 = get_max(x)
end_time <- Sys.time()
t1 <-start_time - end_time

start_time <- Sys.time()
y2 = which.max(x)
end_time <- Sys.time()
t2 <-start_time - end_time
