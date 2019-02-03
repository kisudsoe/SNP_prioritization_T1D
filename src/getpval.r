getpval = function(dist,mu,method='t.test') {
    if(method=='t.test') pval = t.test(dist,mu=mu)$p.value
    else if(method=='z.test') {
        z.score = z.test(dist,var=var(dist),mu=mu)
        pval = pnorm(abs(z.score),lower.tail=F)*2 # two-tailed
    } else if(method=='approxfun') {
        d = density(dist)#; df = approxfun(d)
        pval = try(integrate(approxfun(d),lower=min(d$x),upper=mu)$value)
        #print(paste0('density= ',df(mu)))
        #print(paste0('pval= ',pval))
        if('try-error'%in% class(pval)) pval=0
        if(!is.na(pval)&pval>0.5) pval = 1-pval
    }
    return(pval)
}

suppressMessages(library(ggplot2))
getdistplot = function(dist,mu,title='') {
    # get label y-value
    d = density(dist)
    af = approxfun(d)
    text_y = af(mu)
    if(is.na(text_y)) text_y=0
    
    # get empirical-pval
    pval = try(integrate(approxfun(d),lower=min(d$x),upper=mu)$value)
    if('try-error'%in% class(pval)) pval=0
    if(!is.na(pval)&pval>0.5) pval = 1-pval
    lab = paste0('p=',round(pval,4))
    
    # draw plot
    df = data.frame(Count=d$x,Density=d$y)
    p = ggplot(df,aes(x=Count,y=Density))+theme_bw()+
            geom_line()+
            geom_vline(aes(xintercept=mu),color='red',size=.5)+
            ggtitle(title)+
            geom_text(x=mu,y=text_y,label=lab)
    return(p)
}