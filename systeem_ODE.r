### 10.1 Lotka Volterra model

model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    dR <- r*R*(1 - R/K) - a*R*N
    dN <- c*a*R*N - delta*N
    return(list(c(dR, dN)))
  })
}
p <- c(r=1,K=1,a=1,c=1,delta=0.5) # p is a named vector of parameters
s <- c(R=1,N=0.01) # s is the state

run()
plane()
plane(xmin = -0.001, ymin = -0.001)
plane(xmin = -0.001, ymin = -0.001, portrait = TRUE)
newton(c(R=0.5, N=0.5), plot = TRUE)

### 10.2 The LAC operon
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    R = 1/(1+A^n)
    dA = M*L - delta*A - v*M*A
    dM = c0 + c*(1-R) - d*M
    return(list(c(dA, dM)))
  })
}
p <- c(L=1,c=1,c0=0.05,d=1,delta=0.2,n=5,v=0.25)
s <- c(A=0,M=0)
plane(xmax=4)
low <- newton(s,plot=T)
mid <- newton(c(A=0.8,M=0.2),plot=T)
hig <- newton(c(A=2,M=1),plot=T)
continue(mid,x="L",y="A",xmax=2,ymax=4,add=T)

### 10.1 Tutorial
"1a Dat is het verloop R en N door de tijd met een bepaald startpunt"
"1b Je start op een zeker punt en er wordt bijgehouden wat het systeem doet"
"1c Dit is de phaseplane, hier zie je de nullclines. Op deze lijnen is 1 van de ODE's 0."
"1d Stabiliseert op R is 0.5 en N = 0. Er is niet genoeg carrying capacity om je predator levend te houden."

### 10.2
"2a, nou hier is het model"
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    
    dR <- r*R*(1 - R/K) - a*R*N/(h+R)
    dN <- c*a*R*N/(h+R) - delta*N
    
    return(list(c(dR, dN)))  
  }) 
}  

p <- c(r=1,K=1,h=0.1,a=0.5,c=1,delta=0.4)
s <- c(R=1,N=0.01)

plane(eps=-0.001)
newton(c(R=0.5, N = 0.5), plot = TRUE)
run(200,traject=TRUE)

"2b dan doen we gewoon"
s <- c(R=0.01,N=0)
run()

"2c this with some differenet parameter values"
"Als de steady state unstable is dan gaan het osscileren, noem je dit een limit cycle?"
p <- c(r=1,K=1,h=0.1,a=0.5,c=1,delta=0.4)
s <- c(R=1,N=0.01)
plane(eps=-0.001)
newton(c(R=0.5, N = 0.5), plot = TRUE)
f <- run(200,traject=TRUE)
"2d"
f <- run(200,traject=TRUE)                 # Werkt niet?
f <- newton(f)                             # Werkt niet
f <- newton(c(R=0.5, N = 0.5), plot = TRUE)
continue(f,x="K", y = "R", xmax =2, ymax = 2, add = T)
"2e, nou hier komt ie hoor"
p <- c(r=1,K=0.9,h=0.1,a=0.5,c=1,delta=0.4)
s <- c(R=1,N=0.01)
plane(eps=-0.001)
" De 0.4 bifurcatie is van er kunnen geen predators zijn naar er zijn wel predators."
" De 0.9 bifurcatie is van stabiel naar instabiel"

"2f"
"Als je de K aanpast dan kan je van geen predatoren naar wel predatoren gaan, of van stabiel naar instabiel"

"2g"
"DE PREDATOREN INCREASEN MEER ALS JE DE K VAN DE PREY OMHOOG GOOIT"

"2h"
"Nou we zoeken dus een Lotka Volterra met een stabiel evenwicht en gebruiken dan run() om mooie plotjes te maken"
p <- c(r=1,K=1.4,h=0.1,a=0.5,c=1,delta=0.45)
s <- c(R=0.01,N=0.01)
plane(eps=-0.001)
run(tmax = 1000)
"2i"
"Ja dat kan, bijvoorbeeld door de deathrate omhoog te gooien zodat de predator zich altijd in soort van steady state bevind"


### 10.3 Bacteria and toxin
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    
    dR = s - w*R - a*N*R/(h+R)
    dN = b*N*R/(h+R) - w*N - x*N*Q
    dQ = p*N*Q - d*Q - w*Q
    
    return(list(c(dR, dN, dQ)))  
  }) 
}  

# Parameters taken from Cornejo et al (2009):
p <- c(s=100,a=1e-7,b=1,h=25,w=0.1,x=5e-6,p=4e-10,d=0.1)
s <- c(R=1000,N=1,Q=1)

plane(xmax=6,ymax=4e10,eps=-0.001,portrait=T)
newton(c(R=3,N=1e10,Q=0),plot=T)
f <- run(ymin=1,ymax=1e10,log="y")

"10.3a"
plane(xmax=6,ymax=4e10,eps=-0.001,portrait=T)
"10.3b"
newton(c(R=3,N=1e10,Q=0),plot=T)
"Het is nu instabiel omdat T niet meer steady is, het heeft een soort van vertraging"
"10.3c"
run(tmax = 3000)
"Dit is geen limit cycle, de oscilaties worden kleiner"
"Het is wel unstable, er blijven oscilaties."
"10.3d ik zou de s veranderen"
p <- c(s=10,a=1e-7,b=1,h=25,w=0.1,x=5e-6,p=4e-10,d=0.1)
run(tmax = 3000)
"10.3e"
"Nou dat gaan we eens onderzoeken met het continue command!"
p <- c(s=100,a=1e-7,b=1,h=25,w=0.1,x=5e-6,p=4e-10,d=0.1)
s <- c(R=1000,N=1,Q=1)
plane(xmax=6,ymax=4e10,eps=-0.001,portrait=T)
f <- newton(c(R=3,N=1e10,Q=0),plot=T)
continue(f,x="s", y = "R", xmax =120, add = T)
p <- c(s=70,a=1e-7,b=1,h=25,w=0.1,x=5e-6,p=4e-10,d=0.1)
run(tmax = 3000)
"Dit is geen plotseling overstap"
"10.3f"
"Dit is het nieuwe model"
model <- function(t, state, parms) {
  with(as.list(c(state,parms)), {
    
    dR = s - w*R - a*N*R/(h+R)
    dN = b*N*R/(h+R) - w*N - x*N*Q
    dQ = p0 + p*N*Q - d*Q - w*Q
    
    return(list(c(dR, dN, dQ)))  
  }) 
}  

# Parameters taken from Cornejo et al (2009):
p <- c(s=100,a=1e-7,b=1,h=25,w=0.1,x=5e-6,p=4e-10,d=0.1, p0 = 0.001)
s <- c(R=1000,N=1,Q=1)
run(tmax = 3000)
plane(xmax=6,ymax=4e10,eps=-0.001,portrait=T)
"Ik heb p0 = 0.001 gekozen, dat verandert niet echt iets"
