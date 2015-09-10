library(makehams)
makehams(A=0.00022, B=2.7e-06, c=1.124, d=2, x=20, w=131,
         radix=1e+05, i=0.05)

# Match results shown in Appendix Tables 
# of Actuarial Mathematics for Life Contingent 
# Risks (2nd edition)

# TableD.1
TableD.1 <- createLifeTable()

print(head(TableD.1,5), row.names = F)
print( TableD.1[59:63,], row.names=F)

# TableD.2  -> 5 year pure endowment period
#           -> 1 * force of interest
TableD.2 <- createInsuranceTable(n = 5)

print(head(TableD.2,5), row.names = F)
print( TableD.2[57:61,], row.names=F)
ax <- (1 - TableD.2$"A[x]") / (0.05/(1.05))
print(head(ax,5))
print(ax[57:61])

# TableD.2  -> 10 year endowment period
#           -> 1 * force of interest
TableD.2 <- createInsuranceTable(n = 10)

print(head(TableD.2,5), row.names=F)
print( TableD.2[57:61,], row.names=F)

# TableD.2  -> 20 year endowment period
#           -> 1 * force of interest
TableD.2 <- createInsuranceTable( n = 20)

print(head(TableD.2,5), row.names=F)
print( TableD.2[57:61,], row.names=F)

# TableD.2  -> 5 year endowment period 
#           -> 2 times the force of interest
TableD.2 <- createInsuranceTable( n = 5 , mt = 2)

print(head(TableD.2,5), row.names=F)
print( TableD.2[57:61,], row.names=F)

# TableD.3  -> 5 year endowment period
#           -> 1 * force of interest
TableD.3 <- createInsuranceTable( n = 5, d=0)

print(head(TableD.3,5), row.names=F)
print( TableD.3[57:61,], row.names=F)
ax <- (1 - TableD.3$Ax) / (0.05/(1.05))
print(head(ax,5))
print(ax[57:61])

# TableD.3    -> 10 year endowment period
#             -> 1 * force of interest
TableD.3 <- createInsuranceTable( n = 10, d = 0)

print(head(TableD.3,5), row.names=F)
print(TableD.3[57:61,], row.names=F)

# TableD.3    -> 20 year endowment period
#             -> 1 * force of interest
TableD.3 <- createInsuranceTable( n = 20, d = 0)

print(head(TableD.3, 5), row.names=F)
print(TableD.3[57:61,], row.names=F)

# TableD.3    -> 5 year endowment period
#             -> 2 * force of interest
TableD.3 <- createInsuranceTable( n = 5, d = 0, mt = 2)

print(head(TableD.3, 5), row.names=F)
print(TableD.3[57:61,], row.names=F)

# Table D.4
gl.a(uxt01, function(t,x=gl.g(x)) ifelse(x+t<35,0.1,
  ifelse(x+t<45,0.05,ifelse(x+t<60,0.02,0))))
gl.a(uxt02, function(t,x=gl.g(x)) t^0*0.001)
gl.a(uxt03, function(t,x=gl.g(x))
  ifelse(x+t<60,0,ifelse(x+t<65,0.1,0)))
gl.a(uxt04, function(t,x=gl.g(x), A=gl.g(A), B=gl.g(B),
  c=gl.g(c)) A + B*c^(x+t))
gl.a(uxt, function(t,x=gl.g(x),...) gl.g(uxt01)(t,x) +
  gl.g(uxt02)(t,x) + gl.g(uxt03)(t,x) + gl.g(uxt04)(t,x))
gl.a(radix,1e+06)
tpxij <- function(t,x=gl.g(x),uxt=gl.g(uxt)) {
  integrate(function(s) tpx(s,x)*uxt(s,x), 0, t)$value
}

p = tpx(0:45,20)
st = data.frame(x = 20:65, lx=gl.g(radix)*p, wx=p*gl.g(radix)*
  sapply(0:45, function(k) tpxij(1,20+k,uxt=gl.g(uxt01))),
  ix=p*gl.g(radix)*sapply(0:45, function(k) tpxij(1,20+k, uxt=gl.g(uxt02))), 
  rx=p*gl.g(radix)*sapply(0:45, function(k) tpxij(1,20+k,uxt=gl.g(uxt03))), 
  dx=p*gl.g(radix)*sapply(0:45, function(k) tpxij(1,20+k,uxt=gl.g(uxt04))))

st[41:46,2:6] = st[41:46,2:6]*0.7
st[46,]$wx = 0; st[46,]$dx = 0; st[46,]$ix = 0
st[46,]$rx = st[46,]$lx

TableD.4 <- st
print(head(TableD.4,5), row.names=F)
print(tail(TableD.4,5), row.names=F)
