setwd("C:/Users/jriver58/Dropbox/Postdoc/Share Emilia/JR scel biogeography")


make.tree<-function(edge,edge.length,tips,names=NULL){
	##Rearrange the branches
	for(j in 1:length(edge.length)){
		desc<-which(edge[,1]==edge[j,2])
		if (length(desc)==2){
			temp<-(j+1):length(edge.length)
			temp<-temp[temp!=desc[1] & temp!=desc[2]]
			temp<-c(1:j,desc,temp)
			edge<-edge[temp,]
			edge.length<-edge.length[temp]
		}
	}
	##Renumber the nodes
	ntaxa<-length(tips)
	edge[sort(match(tips,edge[,2])),2]<-1:ntaxa
	root<-edge[1,1]
	edge[edge==root]<-ntaxa+1
	nodes<-edge[edge[,2]<0,2]
	for(j in 1:length(nodes)) edge[which(edge==nodes[j])]<-j+ntaxa+1	

	if (is.null(names)) names<-paste("T",1:ntaxa,sep="")
	else names<-sample(names,length(names))

	tree<-list(edge=edge,tip.label=names,Nnode=length(nodes)+1,
		edge.length=edge.length)
	class(tree)<-"phylo"
	tree
}

tau<-function(taxa,times,lambdas,mus=0,full=FALSE){
	n<-taxa
	if (times[1]<0) stop("Error: times must be positive")
	if (length(times)>1) for (i in 2:length(times)) 
			if(times[i]<times[i-1]) stop("Error: times must be increasing")
	depths<-times-c(0,times[-length(times)])

	if (length(n)>1){
		startn<-c(2,n[-length(n)])
		for (i in 2:length(n)) 
			if(n[i]<n[i-1]) stop("Error: ns must be increasing")
		if (length(n) != length(depths)) 
			stop("Error: ns and depths must be of same length")
		if (is.null(lambdas))lambdas<-log(n/startn)/depths+mus
	}
 
	if (length(depths) != length(mus)) {
		if (length(mus==1)) mus<-rep(mus,length(depths))
		else stop("Must have same length mus and depths")
	}
	

	if (is.null(lambdas)){
		if (length(depths==1)) lambdas<-log(n/2)/depths+mus
		else stop("Error lambda can not be null 
			for a single n with multiple depths")
	}
	if (length(depths) != length(lambdas)) {
		if (length(lambdas==1)) lambdas<-rep(lambdas,length(depths))
		else stop("Must have same length lambdas and depths")
	}

	samelm= lambdas==mus
	eras<-length(depths)

  ##Calculate max taus for each era		
	rs<-lambdas-mus
	expect<-exp(rs*depths)
	taus<-rep(NA,eras)
	nextEs<-rep(0,eras)
	nextTs<-rep(1,eras)
	for (i in eras:1){
		if (samelm[i]){
			starter<-lambdas[i]*depths[i]*(1-nextEs[i])
			denom<-starter+1
			numE<-starter+nextEs[i]
			numT<-nextTs[i]
			numJ<-lambdas[i]*depths[i]*nextTs[i]
		}
		else {
			ender<-lambdas[i]*nextEs[i]-mus[i]
			facto<-(1-nextEs[i])*expect[i]
			denom<-lambdas[i]*facto+ender
			numE<-mus[i]*facto+ender
			numT<-rs[i]^2*expect[i]*nextTs[i]
			numJ<-lambdas[i]*(expect[i]-1)*nextTs[i]
		}
		taus[i]<-numJ/denom
		if (i>1){
			nextEs[i-1]<-numE/denom
			nextTs[i-1]<-numT/denom^2
		}
	}
	if (full) cbind(taus,nextTs,nextEs)
	else taus
}


rbdtree.n<-function(taxa,times,lambdas,mus=0,names=NULL){
	
	n<-taxa
	if (times[1]<0) stop("Error: times must be positive")
	if (length(times)>1) for (i in 2:length(times)) 
			if(times[i]<times[i-1]) stop("Error: times must be increasing")
	depths<-times-c(0,times[-length(times)])

	if (length(n)>1){
		startn<-c(2,n[-length(n)])
		for (i in 2:length(n)) 
			if(n[i]<n[i-1]) stop("Error: ns must be increasing")
		if (length(n) != length(depths)) 
			stop("Error: ns and depths must be of same length")
		if (is.null(lambdas))lambdas<-log(n/startn)/depths+mus
	}
 
	if (length(depths) != length(mus)) {
		if (length(mus==1)) mus<-rep(mus,length(depths))
		else stop("Must have same length mus and depths")
	}
	

	if (is.null(lambdas)){
		if (length(depths==1)) lambdas<-log(n/2)/depths+mus
		else stop("Error lambda can not be null 
			for a single n with multiple depths")
	}
	if (length(depths) != length(lambdas)) {
		if (length(lambdas==1)) lambdas<-rep(lambdas,length(depths))
		else stop("Must have same length lambdas and depths")
	}

	samelm= lambdas==mus
	eras<-length(depths)

##Calculate max taus for each era		
	rs<-lambdas-mus
	expect<-exp(rs*depths)
	taus<-rep(NA,eras)
	nextEs<-rep(0,eras)
	nextTs<-rep(1,eras)
	for (i in eras:1){
		if (samelm[i]){
			starter<-lambdas[i]*depths[i]*(1-nextEs[i])
			denom<-starter+1
			numE<-starter+nextEs[i]
			numT<-nextTs[i]
			numJ<-lambdas[i]*depths[i]*nextTs[i]
		}
		else {
			ender<-lambdas[i]*nextEs[i]-mus[i]
			facto<-(1-nextEs[i])*expect[i]
			denom<-lambdas[i]*facto+ender
			numE<-mus[i]*facto+ender
			numT<-rs[i]^2*expect[i]*nextTs[i]
			numJ<-lambdas[i]*(expect[i]-1)*nextTs[i]
		}
		taus[i]<-numJ/denom
		if (i>1){
			nextEs[i-1]<-numE/denom
			nextTs[i-1]<-numT/denom^2
		}
	}

##Generate random taus
	Enum<-lambdas*nextEs-mus
	Eden<-1-nextEs

	maxtau<-sum(taus)
	if (length(n)==1) raw.taus<-sort(runif(n-2,max=maxtau))
	else {
		other.tau<-maxtau
		raw.taus<-NULL
		for(i in 1:length(n)){
			other.tau<-other.tau-taus[i]
			raw.taus<-c(raw.taus,sort(runif(n[i]-startn[i],
				max=taus[i]))+other.tau)
		}
	}

##Convert taus to times
	times<-rep(NA,n-2)
	which.tau<-length(taus)
	maxtime<-sum(depths)
	for(i in 1:(n-2))
	{
		while(raw.taus[i]>taus[which.tau]){
			raw.taus<-raw.taus-taus[which.tau]
			maxtime<-maxtime-depths[which.tau]
			which.tau<-which.tau-1
		}
		con.tau<-raw.taus[i]/nextTs[which.tau]
		dent<-lambdas[which.tau]*(1-Eden[which.tau]*con.tau)
		if (samelm[which.tau])tempt<-con.tau/dent
		else{
			numt<-lambdas[which.tau]+Enum[which.tau]*con.tau	
			tempt<-log(numt/dent)/rs[which.tau]
		}
		times[i]<-maxtime-tempt
	}

##Combine nodes and build tree
	times<-sort(times)
	cum<-sum(depths)
	edge<-array(c(rep(0,n),1:n),dim<-c(n,2))
	edge.length<-rep(NA,n)
	start.times<-rep(cum,n)
	bottoms<-1:n
	nextn<- -1
	nexttime<-length(times)
	while (length(bottoms)>2){
		combin<-sample(bottoms,2)
		oldb<-match(combin,edge[,2])

		edge[oldb,1]<-nextn
		edge<-rbind(c(0,nextn),edge)

		edge.length[oldb]<-start.times[oldb]-times[nexttime]
		edge.length<-c(0,edge.length)
		
		start.times<-c(times[nexttime],start.times)

		bottoms[bottoms==combin[1]]<-nextn
		bottoms<-bottoms[bottoms!=combin[2]]

		nextn<-nextn-1
		nexttime<-nexttime-1
	}
	lastnodes<-edge[,1]==0
	edge.length[lastnodes]<-start.times[lastnodes]

	make.tree(edge,edge.length,1:n,names)
}

mrbdtree.n<-function(N,taxa,times,lambdas=NULL,mus=0,names=NULL){
	temp<-replicate(N,rbdtree.n(taxa,times,lambdas,mus,names),simplify=FALSE)
	class(temp)<-"multiPhylo"
	temp
}

ltt.null<-function(taxa,times,lambdas,mus=0,intervals=100,method="image",
		cutoffs=c(0.5,0.75,0.95,0.99),colors=NULL,
		x.lab="Time",y.lab="Extant Lineages",show.legend=TRUE,...)
{
	if (times[1]<0) stop("Error: times must be positive")
	if (length(times)>1) for (i in 2:length(times)) 
			if(times[i]<times[i-1]) 
				stop("Error: times must be increasing")
	depths<-times-c(0,times[-length(times)])

	if (length(taxa)>1) if (length(taxa) != length(times))
		stop("Must have taxa legth 1 or same as times")
	
 	if (length(mus==1)) mus<-rep(mus,length(depths))
	else if (length(depths) != length(mus)) stop("Must have same length mus and depths")
	

	if (is.null(lambdas)){
		if (length(depths==1)) lambdas<-log(n/2)/depths+mus
		else stop("Error lambda can not be null 
			for a single n with multiple depths")
	}
	else if (length(lambdas)==1) lambdas<-rep(lambdas,length(depths))
	else if (length(depths) != length(lambdas)) 
		stop("Must have same length lambdas and depths")

	interval_count<-intervals*depths/times[length(times)]
	interval_round<-floor(interval_count)
	interval_round[depths==0]<-1
	interval_count<-interval_count-interval_round
	extra_sort<-sort(interval_count,decreasing=TRUE)
	extra_count<-intervals-sum(interval_round)
	if (extra_count<0) stop ("Error: increase intervals") 
	else if(extra_count>0){
		for (i in 1:extra_count){
			addto<-which(interval_count==extra_sort[i])[1]
			interval_round[addto]<-interval_round[addto]+1
	}}
	times.use<-1:interval_round[1]*depths[1]/(interval_round[1]+1)
	if (length(times)>1) for (i in 2:length(times))
		times.use<-c(times.use,1:interval_round[i]*
		depths[i]/(interval_round[i]+1)+times[i-1])

	count<-1
	lambdas.use<-mus.use<-NULL
	for (i in 1:intervals)
	{
		while(times.use[i]>times[count]) count<-count+1
		lambdas.use<-c(lambdas.use,lambdas[count])
		mus.use<-c(mus.use,mus[count])
	}


	taus<-tau(taxa[length(taxa)],times.use,lambdas.use,mus.use)
	if (length(taxa)==1){ 
		for ( i in 2:intervals) taus[i]<-taus[i]+taus[i-1]
		taus<-taus/taus[intervals]
		taxa.start<-rep(2,intervals)
		taxa.end<-rep(taxa,intervals)
	}
	else { 
		interval_sum<-0
		taxa.start<-taxa.end<-NULL
		for (j in 1:length(taxa)) {
			interval_sum<-c(interval_sum,
				interval_sum[j]+interval_round[j])
			taxa.end<-c(taxa.end,rep(taxa[j],interval_round[j]))
			if (j==1) taxa.start<-rep(2,interval_round[j])
			else taxa.start<-c(taxa.start,rep(taxa[j-1],interval_round[j]))
			if (interval_round[j]>1)
				for ( i in (interval_sum[j]+2):interval_sum[j+1]) 
					taus[i]<-taus[i]+taus[i-1]
			taus[(interval_sum[j]+1):interval_sum[j+1]]<-
				taus[(interval_sum[j]+1):interval_sum[j+1]]/
				taus[interval_sum[j+1]]
		}
	}

	taxa<-taxa[length(taxa)]
	z<-array(0,dim=c(intervals,taxa))
	if (method=="image" | method== "contour") {
		for (x in 1:intervals) {
			for (i in length(cutoffs):1) {
				perc<-(1-cutoffs[i])/2
				bottom<-qbinom(perc,taxa.end[x]-taxa.start[x],taus[x])
				tops<-qbinom(1-perc,taxa.end[x]-taxa.start[x],taus[x])
				z[x,bottom:tops+taxa.start[x]]<-i
			}
		}
		if (is.null(colors)) colors<-heat.colors(length(cutoffs))
		if (method=="image") image(times.use-times[length(times)],1:taxa,z,
			xlab=x.lab,ylab=y.lab,col=colors,
			breaks=0:length(cutoffs)+0.5,...)
		if (show.legend) legend(-times[length(times)]*0.9,0.9*taxa,
			cutoffs,fill=colors)	
	}
}
#ltt.null(c(11,length(plethodon)+1),c(-20,0)+max(plethodon),lam.pleth,mu.pleth)






require(ape)
require(phytools)
require(plyr)
require(treeman)

##read in tree
load("beastLeache.RData")
mult.tree <- beastLeache
scel.tree <- mult.tree[2:1001]
single.tree <- scel.tree[[1]]

contree <- read.nexus("consensus_tree.nxs")



#make null model 
yule.scel <- yule(single.tree)
btimes.scel <- sort(branching.times(single.tree), decreasing=TRUE)
ltt.null(53, btimes.scel[1], yule.scel$lambda, 0)


## smooth out curbs
ltt.null(53, btimes.scel[1], yule.scel$lambda, 0)

temp <- lapply(1:length(scel.tree), function(x) matrix(NA, nrow=53, ncol=2))

for (i in 1:length(scel.tree)){
temp[[i]] <- ltt.plot.coords(scel.tree[[i]])
}

temp2 <- do.call("rbind.fill", lapply(temp, as.data.frame))
smoothScatter(temp2)


plotTree(contree, node.numbers=T)


# nodes to use: 56:79
tip1<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel1",
    edge.length=1.0,
    Nnode=1)
tip2<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel2",
    edge.length=1.0,
    Nnode=1)
tip3<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel3",
    edge.length=1.0,
    Nnode=1)
tip4<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel4",
    edge.length=1.0,
    Nnode=1)
tip5<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel5",
    edge.length=1.0,
    Nnode=1)
tip6<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel6",
    edge.length=1.0,
    Nnode=1)
tip7<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel7",
    edge.length=1.0,
    Nnode=1)
tip8<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel8",
    edge.length=1.0,
    Nnode=1)
tip9<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel9",
    edge.length=1.0,
    Nnode=1)
tip10<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel10",
    edge.length=1.0,
    Nnode=1)
tip11<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel11",
    edge.length=1.0,
    Nnode=1)
tip12<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel12",
    edge.length=1.0,
    Nnode=1)
tip13<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel13",
    edge.length=1.0,
    Nnode=1)
tip14<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel14",
    edge.length=1.0,
    Nnode=1)
tip15<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel15",
    edge.length=1.0,
    Nnode=1)
tip16<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel16",
    edge.length=1.0,
    Nnode=1)
tip17<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel17",
    edge.length=1.0,
    Nnode=1)
tip18<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel18",
    edge.length=1.0,
    Nnode=1)
tip19<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel19",
    edge.length=1.0,
    Nnode=1)
tip20<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel20",
    edge.length=1.0,
    Nnode=1)
tip21<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel21",
    edge.length=1.0,
    Nnode=1)
tip22<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel22",
    edge.length=1.0,
    Nnode=1)
tip23<-list(edge=matrix(c(2,1),1,2),
    tip.label="Scel23",
    edge.length=1.0,
    Nnode=1)
	
class(tip1)<-"phylo"
class(tip2)<-"phylo"
class(tip3)<-"phylo"
class(tip4)<-"phylo"
class(tip5)<-"phylo"
class(tip6)<-"phylo"
class(tip7)<-"phylo"
class(tip8)<-"phylo"
class(tip9)<-"phylo"
class(tip10)<-"phylo"
class(tip11)<-"phylo"
class(tip12)<-"phylo"
class(tip13)<-"phylo"
class(tip14)<-"phylo"
class(tip15)<-"phylo"
class(tip16)<-"phylo"
class(tip17)<-"phylo"
class(tip18)<-"phylo"
class(tip19)<-"phylo"
class(tip20)<-"phylo"
class(tip21)<-"phylo"
class(tip22)<-"phylo"
class(tip23)<-"phylo"
# attach to any node 
contree2<-bind.tree(contree,tip1,where=56)
contree2<-bind.tree(contree2,tip2,where=58)
contree2<-bind.tree(contree2,tip3,where=61)
contree2<-bind.tree(contree2,tip4,where=64)
contree2<-bind.tree(contree2,tip5,where=69)
contree2<-bind.tree(contree2,tip6,where=75)
contree2<-bind.tree(contree2,tip7,where=81)
contree2<-bind.tree(contree2,tip8,where=86)
contree2<-bind.tree(contree2,tip9,where=59)
contree2<-bind.tree(contree2,tip10,where=62)
contree2<-bind.tree(contree2,tip11,where=71)
contree2<-bind.tree(contree2,tip12,where=73)
contree2<-bind.tree(contree2,tip13,where=83)
contree2<-bind.tree(contree2,tip14,where=52)
contree2<-bind.tree(contree2,tip15,where=91)
contree2<-bind.tree(contree2,tip16,where=77)
contree2<-bind.tree(contree2,tip17,where=94)
contree2<-bind.tree(contree2,tip18,where=102)
contree2<-bind.tree(contree2,tip19,where=42)
contree2<-bind.tree(contree2,tip20,where=33)
contree2<-bind.tree(contree2,tip21,where=96)
contree2<-bind.tree(contree2,tip22,where=83)
contree2<-bind.tree(contree2,tip23,where=89)

plotTree(contree2)


	

##start plotting
## first let's make the 2 figs indepedently

## panel 1
ltt.null(53, btimes.scel[1], yule.scel$lambda, 0)
ltt.lines(contree2)

## panel 2
smoothScatter(temp2)
ltt.lines(contree)

#now make them together
par(mfrow=c(1,2))
ltt.null(53, btimes.scel[1], yule.scel$lambda, 0, legend=F)
ltt.lines(contree2)

smoothScatter(temp2)
ltt.lines(contree2)






(3) A value for Pybus & Harvey's "gamma" statistic of -1.2577, p-value = 0.05.
    early burst of speciation, significant






## extra bits we don't need

for (i in 1:length(scel.tree){
ltt.lines(scel.tree[[i]])
}





