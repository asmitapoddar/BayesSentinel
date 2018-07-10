# Noyau Epanechnikov
K_E <- function(u,h){
  (1-(u/h)^2)*(u/h <= 1)
}

# Noyau Gaussien
K_G <- function(u,h){
  exp(-(u^2)/(2*h))
}

# Noyau exponentielle
K_Exp <- function(u,h){
  exp(-(u)/(h))
}

# Noyau Uniforme
K_U <- function(u,h){
  (u/h <= 1)
}

# Noyau Quadratique
K_Q <- function(u,h){
  (1-(u/h)^2)^2*(u/h <= 1)
}

# Noyau Circulaire
K_C <- function(u,h){
  cos(pi/2*u/h)*(u/h <= 1)
}

# Noyau Triangulaire
K_T<- function(u,h){
  (1-u/h)*(u/h <= 1)
}

#Rational Quadratique
K_RQ <- function(u,h){
  1 - (u^2)/((u^2)+h)
}

#Inverse Multiquadratic
K_IMQ <- function(u,h){
  1/sqrt((u)^2 + h^2)
}
