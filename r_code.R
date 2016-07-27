#read in the environmental covariates
n = 16 #num covariates (inc. interactions etc.)
x1 = as.matrix(read.table("data/x1.asc")) #Annual Mean Temperature
x2 = as.matrix(read.table("data/x2.asc")) #Annual Precipitation
x3 = as.matrix(read.table("data/x3.asc")) #Temperature Seasonality
x4 = as.matrix(read.table("data/x4.asc")) #Precipitation of Wettest Quarter
x5 = as.matrix(read.table("data/x5.asc")) #Precipitation of Driest Quarter
x6 = as.matrix(read.table("data/x6.asc")) #Precipitation of Warmest Quarter
x7 = as.matrix(read.table("data/x7.asc")) #Precipitation of Coldest Quarter
z1 = as.matrix(read.table("data/z1.asc")) #Human Footprint
z2 = as.matrix(read.table("data/z2.asc")) #iNat Footprint

#read in the state map
states = as.matrix(read.table("data/states.asc")) #States

#read in the species data
obs = read.csv("data/obs.csv")
lt = t(read.csv("data/listings.csv")); dim(lt) = NULL

#Set up the env
x = t(read.csv("data/x.csv")); dim(x) = NULL
y = t(read.csv("data/y.csv")); dim(y) = NULL
by = 0.1666666666666666

mask = matrix(0,length(x),length(y))
for(i in lt){
  mask[which(states==i)] = 1
}
dim(mask) = NULL
maskdat = cbind(expand.grid(x,y),mask)
lims = apply(maskdat[which(maskdat[,3]==1),1:2],2,range)
lims = lims+matrix(c(-5,5,-5,5),2,2)
cols = c("#000000","#070400","#0E0800","#160C00","#1D1000","#241400","#2C1800","#331C00","#3B2100","#422500","#492900","#512D00","#583100","#603500","#673900","#6E3D00","#764200","#7D4600","#844A00","#8C4E00","#935200","#9B5600","#A25A00","#A95E00","#B16300","#B86700","#C06B00","#C76F00","#CE7300","#D67700","#DD7B00","#E58000","#E58408","#E68810","#E78C18","#E89020","#E99429","#EA9831","#EA9C39","#EBA041","#ECA44A","#EDA852","#EEAD5A","#EFB162","#EFB56A","#F0B973","#F1BD7B","#F2C183","#F3C58B","#F4C994","#F4CD9C","#F5D1A4","#F6D6AC","#F7DAB4","#F8DEBD","#F9E2C5","#F9E6CD","#FAEAD5","#FBEEDE","#FCF2E6","#FDF6EE","#FEFAF6","#FFFFFF")

#run the model
pomat = matrix(NA,0,n)
for(i in 1:dim(obs)[1]){
  slot_y = intersect(which( y < obs$latitude[i] ), which( (y + by) > obs$latitude[i] ))
  slot_x = intersect(which( x < obs$longitude[i] ), which( (x + by) > obs$longitude[i] ))
  pomat = rbind(pomat,c(x1[slot_x,slot_y], x2[slot_x,slot_y], x3[slot_x,slot_y], x4[slot_x,slot_y], x5[slot_x,slot_y], x6[slot_x,slot_y], x7[slot_x,slot_y], x1[slot_x,slot_y]^2, x2[slot_x,slot_y]^2, x3[slot_x,slot_y]^2, x4[slot_x,slot_y]^2, x5[slot_x,slot_y]^2, x6[slot_x,slot_y]^2, x7[slot_x,slot_y]^2, z1[slot_x,slot_y], z2[slot_x,slot_y]))
}
pos = as.data.frame(pomat[apply(pomat,1,function(x){!any(is.na(x))}),])
colnames(pos) = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "z1", "z2")

crit = intersect(intersect(intersect(which(!is.na(x1)), which(!is.na(x2))), which(!is.na(z1))),which(mask==1))
stack = cbind(x1[crit], x2[crit], x3[crit], x4[crit], x5[crit], x6[crit], x7[crit], x1[crit]^2, x2[crit]^2, x3[crit]^2, x4[crit]^2, x5[crit]^2, x6[crit]^2, x7[crit]^2, z1[crit], z2[crit])
p = 1000 #background sample
if(dim(stack)[1] < p){
  p = dim(stack[1])
}
 
p1 = cbind(rep(1,times=dim(pos)[1]),pos)
colnames(p1) = c("y", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "z1", "z2")
stack3 = stack[sample(1:dim(stack)[1],p),]
p2 = cbind(rep(0,times=dim(stack3)[1]),stack3)
colnames(p2) = c("y", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14", "z1", "z2")
pdats = rbind(p1,p2)

pmod = glm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+z1+z2, data=pdats, family=poisson(link = "log"))

#create a distribution map
coefs = summary(pmod)$coefficients[,1]
pred = exp(stack[,c(1:14)]%*%coefs[c(2:15)])
predmat = matrix(x1*0,length(x),length(y))
predmat[crit] = pred
png("map.png")
image(x,y,predmat, asp=1, xlim=lims[,1], ylim=lims[,2], xlab="", ylab="", main="", col=cols)
dev.off()
