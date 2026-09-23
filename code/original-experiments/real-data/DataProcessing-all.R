
load("37070-0001-Data.rda")

dt <- da37070.0001

table(dt$SCHID)

raw.list <- lcc.list <- list()


for(MM in unique(dt$SCHID)){

df <- dt[dt$SCHID==MM,]
dim(df)


df$SCHID <- as.character(df$SCHID)
df$SCHIDW2 <- as.character(df$SCHIDW2)

df$ID <- as.character(df$ID)
id.tab <- table(df$ID)

df$ID <- gsub(" ","",df$ID)

df$IDN <- as.character(df$IDN)
df$IDN <- gsub(" ","",df$IDN)
df$ID <- as.numeric(df$ID)
df$IDN <- as.numeric(df$IDN)


df <- df[!is.na(df$ID),]

df <- df[as.numeric(df$ID)>0,]
df <- df[as.numeric(df$ID)<999,]

ID.tab <- table(df$ID)
ID.tab

final.ID <- as.numeric(names(ID.tab)[which(ID.tab==1)])
n <- max(final.ID)
sub.index <- unlist(lapply(df$ID,function(x){any(final.ID==x)}))
df <- df[sub.index,]

rownames(df) <- as.character(df$ID)

B1 <- DA1 <- A1 <- A <- DA <- B <- matrix(0,n,n)
rownames(B1) <- colnames(B1) <- rownames(B) <- colnames(B) <- rownames(A1) <- colnames(A1) <- rownames(A) <- colnames(A) <- rownames(DA) <- colnames(DA) <- as.character(1:n) 
for(k in 1:nrow(df)){
  i <- df$ID[k]
  STs <- unique(df[k,c(234:243)])
  valid.ST <- which(unlist(lapply(STs,function(x){any(final.ID==x)})))
  STs <- as.numeric(STs[valid.ST])
  A[i,STs] <- 1
  DA[i,STs] <- 1
  A[STs,i] <- 1
  STW1s <- unique(df[k,c(66:75)])
  valid.STW1 <- which(unlist(lapply(STW1s,function(x){any(final.ID==x)})))
  STW1s <- as.numeric(STW1s[valid.STW1])
  A1[i,STW1s] <- 1
  DA1[i,STW1s] <- 1
  A1[STW1s,i] <- 1
  BFs <- unique(df[k,c(249,250)])
  valid.BF <- which(unlist(lapply(BFs ,function(x){any(final.ID==x)})))
  BFs <- as.numeric(BFs[valid.BF])
  B[i,BFs] <- 1
  BFs <- unique(df[k,c(81,82)])
  valid.BF <- which(unlist(lapply(BFs ,function(x){any(final.ID==x)})))
  BFs <- as.numeric(BFs[valid.BF])
  B1[i,BFs] <- 1
}

diag(A) <- 0
diag(DA) <- 0
diag(B) <- 0
diag(A1) <- 0
diag(B1) <- 0
norm(A-t(A),"F")
norm(A-t(DA),"F")

summary(rowSums(A))
summary(rowSums(DA))
summary(rowSums(B))

A2 <- A
B2 <- B



B1.sym <- B1+t(B1)
B1.sym[B1.sym>0] <- 1

B2.sym <- B2+t(B)
B2.sym[B2.sym>0] <- 1

AB1 <- A1+B1.sym
AB1[AB1>0] <- 1

AB2 <- A2+B2.sym
AB2[AB2>0] <- 1


A <- (AB1+AB2)/2


raw.list[[MM]] <- list(df=df,A2=A2,A1=A1,DA=DA,B=B2,B1=B1,DA1=DA1,A=A)

library(igraph)

g <- graph.adjacency(A,mode="undirected")

cl <- clusters(graph = g,mode = "weak")
cl$csize

lcc.index <- which(cl$membership==which.max(cl$csize))


A <- A[lcc.index,lcc.index]
DA <- DA[lcc.index,lcc.index]
B2 <- B2[lcc.index,lcc.index]
B1 <- B1[lcc.index,lcc.index]
A1 <- A1[lcc.index,lcc.index]
A2 <- A2[lcc.index,lcc.index]

df <- df[rownames(A),]

dim(A)
dim(df)

g <- graph.adjacency(A,mode="undirected")


plot(g,vertex.label=NA,vertex.size=4)

lcc.list[[MM]] <- list(df=df,A2=A2,A1=A1,DA=DA,B=B2,B1=B1,DA1=DA1,A=A)


}

save(raw.list,lcc.list,file="ProcessedData.Rda")


for(k in 1:60){
  DA <- lcc.list[[k]]$DA
  B <- lcc.list[[k]]$B
  
  write.table(DA,file=paste("School_",k,"_SpendingTimeAdj.txt",sep=" "))
  write.table(B,file=paste("School_",k,"_BestFriendAdj.txt",sep=" "))
}
