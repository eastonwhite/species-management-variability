#--Distrav--
# Calculates the distance travelled in a landscape.
# N vector of abundances
# th threshold abundance for detection
#
Distrav <- function(N,th){
    max(seq(1,length(N)) * (N >= th)) - 1
}
