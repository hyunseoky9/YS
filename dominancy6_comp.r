genotype2index <- function(genotype, new_al_num) {
# converts the genotype to the index in new population array 
  if(min(genotype)==1){
    index = max(genotype)
  } else {
    index = sum(seq(new_al_num,1)[1:(min(genotype)-1)]) +
    abs(genotype[2]-genotype[1]) + 1
  }
  return(index)
}