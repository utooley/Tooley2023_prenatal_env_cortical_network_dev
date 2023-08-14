library(tidyr)
library(mgcv)
library(gratia)
library(tidyverse)
library(dplyr)
library(lme4)
library(pbkrtest)

#### RESIDUAL PLOT FOR GAMM INTERACTION, WITH SPAGHETTI LINES ####
# Reworked from code from Dr. Bart Larsen
resid_plot_int <- function(modobj,term,int.term=NULL,add.intercept=FALSE, group.var="modid", low_color="#FAA820", high_color="#4D8AC8"){
  library(gratia)
  font_size=16
  theme_set(theme_classic(base_family = "sans",base_size = font_size))
  terms=c(term,int.term)
  model_terms <- term_names(modobj)
  # Check the model
  if (any(class(modobj)=="gam")) {
    modobj <- modobj
  } else if (class(modobj$gam)=="gam") {
    modobj <- modobj$gam
  } else {
    stop("Can't find a gam object")
  }
  
  df <- modobj$model%>%select(model_terms)
  mod.intercept <- modobj$coefficients["(Intercept)"]
  
  variables <- attr(modobj$terms,"term.labels")
  variableClasses <- attr(modobj$terms,"dataClasses")
  thisPred <- list()
  
  for (v in c(1:length(variables))) {
    thisVar <- variables[[v]]
    thisClass <- variableClasses[thisVar]
    if (thisVar == term) {
      thisPred[[term]] = df[,term]
    } else if (thisVar == int.term) {
      switch ( thisClass,
               "numeric" = {thisPred[[thisVar]] = quantile(df[,thisVar],probs=c(.10,.90))},
               "factor" = {thisPred[[thisVar]] = levels(df[,thisVar])},
               "ordered" = {thisPred[[thisVar]] = levels(df[,thisVar])}
      )
    } else {
      switch (thisClass,
              "numeric" = {thisPred[[thisVar]] = mean(df[,thisVar])},
              "factor" = {thisPred[[thisVar]] = levels(df[,thisVar])[[1]]},
              "ordered" = {thisPred[[thisVar]] = levels(df[,thisVar])[[1]]}
      )
    }
  }
  
  variable_expressions <- list()
  for (variable_name in variables) {
    #print(variable_name)
    variable_expression_str <- sprintf("%s = %s", variable_name, deparse(substitute(thisPred[[variable]], list(variable=variable_name))))
    #print(variable_expression_str)
    variable_expression <- eval(parse(text = variable_expression_str))
    variable_expressions[[variable_name]] <- variable_expression
  }
  
  #variable_expressions
  newdf <- do.call(expand.grid,variable_expressions)
  
  pterms <- predict(modobj,type = "terms",se.fit = TRUE,newdata = newdf) #this is trying to get the prediction line
  pterms.full <- predict(modobj,type = "terms",se.fit = TRUE) #full is the real data
  
  pterms.fit <- pterms$fit
  pterms.sefit <- pterms$se.fit
  pterms.full.fit <- pterms.full$fit
  pterms.full.sefit <- pterms.full$se.fit
  
  pterms.df <- data.frame(pterms.fit)%>%select(matches(terms))%>%
    mutate(fit=rowSums(across(where(is.numeric))))%>% #sum across the fit for main term and interaction term
    select(fit)%>%
    # rename(fit:=!!term)%>%
    # cbind(data.frame(pterms.sefit)%>%select(matches(term))%>%rename(se.fit:=!!term))%>%
    cbind(data.frame(pterms.sefit)%>%select(matches(terms))%>%mutate(se.fit=rowSums(across(where(is.numeric)))))%>%
    mutate(upr = fit + 1.96*se.fit,
           lwr = fit - 1.96*se.fit)%>%
    #mutate(partial.residuals=fit+resid(modobj))%>%
    select(fit,se.fit,lwr,upr)
  pterms.df <- cbind(pterms.df,newdf%>%select(matches(terms))%>%rename_all(.funs = function(x) paste0(x,".raw")))
  
  pterms.full.df <- data.frame(pterms.full.fit)%>%select(matches(terms))%>%
    mutate(fit=rowSums(across(where(is.numeric))))%>% #sum across the fit for main term and interaction term
    select(fit)%>%
    # rename(fit:=!!term)%>%
    # cbind(data.frame(pterms.sefit)%>%select(matches(term))%>%rename(se.fit:=!!term))%>%
    cbind(data.frame(pterms.full.sefit)%>%select(matches(terms))%>%mutate(se.fit=rowSums(across(where(is.numeric)))))%>%
    mutate(upr = fit + 1.96*se.fit,
           lwr = fit - 1.96*se.fit)%>%
    mutate(partial.residuals=fit+resid(modobj))%>%
    select(fit,se.fit,lwr,upr, partial.residuals)
  
  partial.residuals.df <- cbind(pterms.full.df,df%>%select(matches(terms))%>%rename_all(.funs = function(x) paste0(x,".raw")))
  group <- modobj$model[,group.var]
  partial.residuals.df <- cbind(partial.residuals.df,group)
  
  if (add.intercept==TRUE) {
    partial.residuals.df$partial.residuals <- partial.residuals.df$partial.residuals+mod.intercept
    pterms.df$fit <-  pterms.df$fit+mod.intercept
    pterms.df$lwr <-  pterms.df$lwr+mod.intercept
    pterms.df$upr <-  pterms.df$upr+mod.intercept
  } else{
  }
  
  ggplot(data=NULL)+
    geom_point(data = partial.residuals.df, aes_string(x=paste0(term,".raw"), y="partial.residuals",
                                                       fill=paste0(int.term,".raw")),
               alpha=.75,color="white",pch=21)+
    geom_line(data = partial.residuals.df, aes_string(x=paste0(term,".raw"), y="partial.residuals", group = "group", color=paste0(int.term,".raw")),alpha = .3) +
    scale_fill_gradient(low = low_color, high = high_color)+
    scale_color_gradient(low = low_color, high = high_color)+
    geom_line(data=filter(pterms.df, names(get(paste0(int.term,".raw")))=="90%" ),aes_string(x= paste0(term,".raw"),y="fit"),
              linewidth=1.5, color=high_color) +
    geom_line(data=filter(pterms.df, names(get(paste0(int.term,".raw")))=="10%" ),aes_string(x= paste0(term,".raw"),y="fit"), color=low_color,
              linewidth=1.5)+
    geom_ribbon(data =filter(pterms.df, names(get(paste0(int.term,".raw")))=="10%"), aes_string(x=paste0(term,".raw"),y="fit",ymin="lwr",ymax="upr"),fill = low_color,
                alpha=.55,color=NA)+
    geom_ribbon(data =filter(pterms.df, names(get(paste0(int.term,".raw")))=="90%"), aes_string(x=paste0(term,".raw"),y="fit",ymin="lwr",ymax="upr"),
                fill = high_color,
                alpha=.55,color=NA)+
    xlab(term)+ylab("partial.residuals")
}

#### RESIDUAL PLOT FOR GAMMS ####
library(scales)
font_size=14
theme_set(theme_classic(base_family = "sans",base_size = font_size))
resid_plot <- function(modobj,term,group.var,int.term=NULL,add.intercept=FALSE, high_color="#4D8AC8", low_color="#FAA820"){
  # Check the model
  if (any(class(modobj)=="gam")) {
    modobj <- modobj
  } else if (class(modobj$gam)=="gam") {
    modobj <- modobj$gam
  } else {
    stop("Can't find a gam object")
  }
  
  terms=c(term,int.term)
  df <- modobj$model
  mod.intercept <- modobj$coefficients["(Intercept)"]
  pterms <- predict(modobj,type = "terms",se.fit = TRUE)
  
  if (add.intercept==TRUE) {
    pterms.fit <- pterms$fit+mod.intercept
  } else{
    pterms.fit <- pterms$fit
  }
  pterms.sefit <- pterms$se.fit
  
  colnames(pterms.fit) <- gsub(x = colnames(pterms.fit),pattern = "s\\(",replacement = "")%>%
    gsub(pattern = "\\)",replacement = "")
  colnames(pterms.sefit) <- gsub(x = colnames(pterms.sefit),pattern = "s\\(",replacement = "")%>%
    gsub(pattern = "\\)",replacement = "")
  
  pterms.df <- data.frame(pterms.fit)%>%select(matches(terms))%>%
    mutate(fit=rowSums(across(where(is.numeric))))%>% #sum across the fit for main term and interaction term
    select(fit)%>%
    # rename(fit:=!!term)%>%
    # cbind(data.frame(pterms.sefit)%>%select(matches(term))%>%rename(se.fit:=!!term))%>%
    cbind(data.frame(pterms.sefit)%>%select(matches(terms))%>%mutate(se.fit=rowSums(across(where(is.numeric)))))%>%
    mutate(upr = fit + 1.96*se.fit,
           lwr = fit - 1.96*se.fit)%>%
    mutate(partial.residuals=fit+resid(modobj))%>%
    select(fit,se.fit,lwr,upr,partial.residuals)
  
  partial.residuals.df <- cbind(pterms.df,df%>%select(matches(terms))%>%rename_all(.funs = function(x) paste0(x,".raw")))
  group <- df[,group.var]
  partial.residuals.df <- cbind(partial.residuals.df,group)
  
  # partial.residuals.df<-data.frame(pterms.fit)%>%
  #   mutate(across(.cols = everything(),.fns = function(x){x+resid(modobj)}))%>%
  #   cbind(df%>%select(matches(terms))%>%rename_all(.funs = function(x) paste0(x,".raw")))
  
  # plot.df <- cbind(partial.residuals,pterms.df)
  
  #make a newdata dataframe of values to predict over for the fit and CI lines, rather than the real data
  ## Generate custom line plot
  np <- 10000 #number of predicted values
  df <- modobj$model
  
  theseVars <- attr(modobj$terms,"term.labels")
  varClasses <- attr(modobj$terms,"dataClasses")
  thisResp <- as.character(modobj$terms[[2]])
  
  if (!is.null(int.term)) {
    # We will produce and interaction plot
    if (is.numeric(partial.residuals.df[,paste0(int.term,".raw")])){
      partial.residuals.df$int.split <- as.factor(ifelse(partial.residuals.df[,paste0(int.term,".raw")] < median(partial.residuals.df[,paste0(int.term,".raw")]), 0, 1))
    }
    if (!any(grepl(x=as.character(modobj$formula),pattern = int.term))) {
      warning("int.term not recognized in modobj formula!")
      return()
    }
    switch (varClasses[int.term],
            "numeric" = {
              q <- quantile(df[,int.term],probs = c(.05,.95)) #pick 10% and 90% to plot
              bigq <- q[[2]]
              smallq <- q[[1]]
              values <- c(bigq,smallq)
              labs <- c(sprintf("high (%1.2f)",bigq),sprintf("low (%1.2f)",smallq))
              
              q <-quantile(rescale(df[,int.term],c(0,1)),probs = c(0,.5,1))
              limit_values <- c(q[[1]],q[[length(q)]])
              midpoint_val <- unname(q[[2]])
              cbar_vals <- unname(q)
              
              theseLabs = rep(values,each = np)
              grad_fill = T
            },
            "factor" = {
              labs <- levels(df[,int.term])
              values <- levels(df[,int.term])
              theseLabs = rep(values,each = np)
              grad_fill = F
            },
            "ordered" = {
              labs <- levels(df[,int.term])
              values <- levels(df[,int.term])
              theseLabs = ordered(rep(values,each = np),levels = values)
              grad_fill = F
            }
    )
    
    labPred <- data.frame(init = rep(0,np*length(labs)))
    labPred[,int.term] = theseLabs
    labPred$lab = rep(labs,each = np)
    labPred <- labPred[,names(labPred) !="init"]
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == term) {
        thisPred[,term] = seq(min(df[,term],na.rm = T),max(df[,term],na.rm = T), length.out = np)
      } else if (v == int.term) {
        next
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    
    thisPred <- thisPred %>% select(-init)
    thisPred <- do.call("rbind", replicate(length(labs), thisPred, simplify = FALSE))
    
    pred <- cbind(labPred,thisPred)
    p<-data.frame(predict(modobj,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    pred[,group.var] = NA #these columns have to exist in the dataframe for plotting
    pred[,thisResp] = 1 #these columns have to exist in the dataframe for plotting
    
    #shift the lines up to meet those of the partial residuals
    difference= mean(partial.residuals.df$partial.residuals) - mean(pred$fit)
    pred$fit <- pred$fit+difference
    pred$selo <- pred$selo+difference
    pred$sehi <- pred$sehi+difference
  }
  #leave this down here
  if (is.null(int.term)) {
    ggplot(partial.residuals.df,aes_string(x=paste0(term,".raw"),y="partial.residuals"))+
      geom_point(alpha=.75,color="white",pch=21,fill="gray")+
      geom_line(aes_string(group = "group"),alpha = .3) +
      geom_ribbon(aes(y=fit,ymin=lwr,ymax=upr),fill="black",alpha=.25,color=NA)+
      geom_line(aes(y=fit),color="black",linewidth=1.5)+
      xlab(term)+ylab("partial.residuals")
  } else{
    ggplot()+
      geom_point(data=partial.residuals.df,aes_string(x=paste0(term,".raw"),
                                                      y="partial.residuals",
                                                      color="int.split",
                                                      fill="int.split"), alpha=.75,color="white",pch=21)+
      geom_line(data=partial.residuals.df,aes_string(x=paste0(term,".raw"),
                                                     y="partial.residuals",
                                                     color="int.split",
                                                     group = "group"),alpha = .3) +
      #geom_ribbon(aes(y=fit,ymin=lwr,ymax=upr),alpha=.55,color=NA)+
      #geom_line(aes(y=fit),linewidth=1.5)+
      geom_ribbon(data = pred, aes_string(x = term , ymin = "selo",ymax = "sehi", fill = "lab"),alpha=.55,color=NA) +
      geom_line(data = pred,aes_string(x = term, y = "fit",col = "lab"), linewidth=0.7) +
      scale_fill_manual(values = c(low_color, high_color,high_color,low_color))+
      scale_color_manual(values = c(low_color, high_color,high_color, low_color))+
      #scale_color_brewer(type = "qual",palette = "Set2")+
      #scale_fill_brewer(type = "qual",palette = "Set2")+
      xlab(term)+ylab("partial.residuals")
  }
}

#### PLOTTING METHOD FOR DERIVATIVES OF GAMS ####
## Used with permission from Dr. Bart Larsen
library(broom)
library(kableExtra)
library(RColorBrewer)
font_size <- 12
#theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2
visualize_model <- function(modobj,smooth_var, int_var = NULL ,group_var, plabels = NULL,check_diagnostics = F,derivative_plot = F, resid_spaghetti=T){
  this_font_size = font_size*1.25
  theme_set(theme_classic(base_family = "sans",base_size = font_size))
  if (any(class(modobj)=="gam")) {
    model <- modobj
  } else if (class(modobj$gam)=="gam") {
    model <- modobj$gam
  } else {
    stop("Can't find a gam object to plot")
  }
  s<-summary(model)
  
  ## Generate custom line plot
  np <- 10000 #number of predicted values
  df = model$model
  
  theseVars <- attr(model$terms,"term.labels")
  varClasses <- attr(model$terms,"dataClasses")
  thisResp <- as.character(model$terms[[2]])
  
  if (!is.null(int_var)) {
    # We will produce and interaction plot
    if (!any(grepl(x=as.character(model$formula),pattern = int_var))) {
      warning("int_var not recognized in model formula!")
      return()
    }
    switch (varClasses[int_var],
            "numeric" = {
              q <- quantile(df[,int_var],probs = c(.05,.95)) #pick 10% and 90% to plot
              bigq <- q[[2]]
              smallq <- q[[1]]
              values <- c(bigq,smallq)
              labs <- c(sprintf("high (%1.2f)",bigq),sprintf("low (%1.2f)",smallq))
              
              q <-quantile(scales::rescale(df[,int_var],c(0,1)),probs = c(0,.5,1))
              limit_values <- c(q[[1]],q[[length(q)]])
              midpoint_val <- unname(q[[2]])
              cbar_vals <- unname(q)
              
              theseLabs = rep(values,each = np)
              grad_fill = T
            },
            "factor" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = rep(values,each = np)
              grad_fill = F
            },
            "ordered" = {
              labs <- levels(df[,int_var])
              values <- levels(df[,int_var])
              theseLabs = ordered(rep(values,each = np),levels = values)
              grad_fill = F
            }
    )
    
    labPred <- data.frame(init = rep(0,np*length(labs)))
    labPred[,int_var] = theseLabs
    labPred$lab = rep(labs,each = np)
    labPred <- labPred[,names(labPred) !="init"]
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else if (v == int_var) {
        next
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    
    thisPred <- thisPred %>% select(-init)
    thisPred <- do.call("rbind", replicate(length(labs), thisPred, simplify = FALSE))
    
    pred <- cbind(labPred,thisPred)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    pred[,group_var] = NA #these columns have to exist in the dataframe for plotting
    pred[,thisResp] = 1 #these columns have to exist in the dataframe for plotting
    
    df$realp <- as.vector(predict(model,df)) #added this here
    
    low_color = "#FAA820" #6EA802 green
    high_color = "#4D8AC8" #C528F4 purple
    high_line = "#FAA820"
    low_line = "#4D8AC8"
    
    if (grad_fill == T) {
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var)) +
        geom_point(alpha = 0.55,stroke = 0, size = point_size) + geom_line(aes_string(group = group_var),alpha = .5) +
        scale_color_gradientn(colors = c(low_color,high_color), values = cbar_vals,name = "") +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = "lab"),alpha = .18, linetype = 0) +
        scale_fill_manual(values = c(high_color,low_color)) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",group = "lab"),size = line_size) +
        labs(title = plabels)+
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp, color = int_var)) +
        geom_point(alpha = .35,stroke = 0, size = point_size) + geom_line(aes_string(group = group_var),alpha = .3) +
        scale_color_brewer(type = "qual",palette = "Set1",direction = 1) +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi", fill = int_var),alpha = .5, linetype = 0) +
        scale_fill_brewer(type = "qual",palette = "Set1",direction = 1) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit",color = int_var),size = line_size) +
        labs(title = plabels)+
        theme(plot.title = element_text(hjust = 0.5))
    }
  } else {
    int_var=""
    # No interaction variable, just produce a single line plot
    thisPred <- data.frame(init = rep(0,np))
    
    for (v in c(1:length(theseVars))) {
      thisVar <- theseVars[[v]]
      thisClass <- varClasses[thisVar]
      if (thisVar == smooth_var) {
        thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np)
      } else {
        switch (thisClass,
                "numeric" = {thisPred[,thisVar] = median(df[,thisVar])},
                "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]},
                "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}
        )
      }
    }
    pred <- thisPred %>% select(-init)
    p<-data.frame(predict(model,pred,se.fit = T))
    pred <- cbind(pred,p)
    pred$selo <- pred$fit - 2*pred$se.fit
    pred$sehi <- pred$fit + 2*pred$se.fit
    pred[,group_var] = NA
    pred[,thisResp] = 1
    if (resid_spaghetti==F) {
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp)) +
        geom_point(alpha = .3,stroke = 0, size = point_size) + geom_line(aes_string(group = group_var),alpha = .3) +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi"),alpha = .5, linetype = 0) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit"),size = line_size) +
        labs(title = plabels)
    } else {
      #still need to figure out how to make this predict the values for individual participants timepoints!
      p1 <- ggplot(data = df, aes_string(x = smooth_var,y = thisResp)) +
        geom_ribbon(data = pred,aes_string(x = smooth_var , ymin = "selo",ymax = "sehi"),alpha = .5, linetype = 0) +
        geom_line(data = pred,aes_string(x = smooth_var, y = "fit"),size = line_size) +
        labs(title = plabels)
    }
    
  }
  
  if (derivative_plot == T) {
    # We will add a bar that shows where the derivative is significant.
    # First make some adjustments to the line plot.
    p1<- p1+theme(text = element_text(size=this_font_size),
                  axis.text = element_text(size = this_font_size),
                  axis.title.y = element_text(size = this_font_size),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  legend.text = element_text(size = this_font_size),
                  legend.title = element_text(size = this_font_size),
                  axis.title = element_text(size = this_font_size),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(fill = "transparent",colour = NA),
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  plot.margin = unit(c(.2, .2, 0, .2), "cm")) #Top, left,Bottom, right
    scatter <- list(p1)
    
    # Now add the plots using the derivative plotting function
    if (any(grepl(x = row.names(s$s.table),pattern =  ":") & grepl(x=row.names(s$s.table),pattern = int_var))) {
      # Factor levels separately if there is an interaction in the model.
      f<-formula(model) # current formula
      fterms <- terms(f)
      fac <- attr(fterms, "factors")
      idx <- which(as.logical(colSums(fac[grep(x=row.names(fac),pattern = int_var),])))
      new_terms <- drop.terms(fterms, dropx = idx, keep.response = TRUE)
      new_formula <- formula(new_terms) # Formula without any interaction terms in the model.
      
      #add derivative gradients for each level of the factor
      num_levels <- length(levels(df[,int_var]))
      level_colors <- suppressWarnings(RColorBrewer::brewer.pal(num_levels,"Set1")) #use the same palette as the line plot
      plotlist = vector(mode = "list",length = num_levels+1) # we will be building a list of plots
      plotlist[1] = scatter # first the scatter plot
      
      for (fcount in 1:num_levels) {
        this_level <- levels(df[,int_var])[fcount]
        df$subset <- df[,int_var] == this_level
        df$group_var <- df[,group_var]
        this_mod <- gamm(formula = new_formula,data = df,subset = subset,random=list(group_var=~1))
        # this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        this_d <- get_derivs_and_plot(modobj = this_mod,smooth_var = smooth_var,low_color = "white",hi_color = level_colors[fcount])
        
        if (fcount != num_levels & fcount != 1){
          # get rid of redundant junk
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
          this_d$theme$legend.text = element_blank()
        }
        if (fcount == 1) {
          this_d$theme$axis.title = element_blank()
          this_d$theme$axis.text.x = element_blank()
          this_d$theme$axis.ticks=element_blank()
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.text = element_blank()
        }
        if (fcount == num_levels) {
          this_d$theme$legend.background=element_blank()
          this_d$theme$legend.box.background = element_blank()
          this_d$theme$legend.key = element_blank()
          this_d$theme$legend.title = element_blank()
        }
        this_d$labels$fill=NULL
        plotlist[fcount+1] = list(this_d)
      }
      pg<-plot_grid(rel_heights = c(16*num_levels,rep(num_levels,num_levels-1),3*num_levels),plotlist = plotlist,align = "v",axis = "lr",ncol = 1)
      final_plot <- pg
      print(final_plot)
    } else {
      # No need to split
      d1 <- get_derivs_and_plot(modobj = modobj,smooth_var = smooth_var)
      scatter <- list(p1)
      bar <- list(d1)
      allplots <- c(scatter,bar)
      pg<-plot_grid(rel_heights = c(16,3),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
      final_plot <- pg
      print(final_plot)
    }
    
  }    else {
    # No derivative plot
    p1<- p1+theme(text = element_text(size=font_size),
                  axis.text = element_text(size = font_size),
                  legend.text = element_text(size = font_size),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  plot.background = element_blank())
    final_plot <- p1
    print(final_plot)
  }
  
  if (check_diagnostics == T) {
    cp <- check(b,
                a.qq = list(method = "tnorm",
                            a.cipoly = list(fill = "light blue")),
                a.respoi = list(size = 0.5),
                a.hist = list(bins = 10))
    print(cp)
  }
  return(final_plot)
}

### function to extract derivative, confidence interval, significance, and plot ###
### This works for signle smooth terms and factor-smooth interactions. Does not work for bivariate smooths.
### This part is still under development.
### If you want to plot derivatives for a by-variable factor model, and you want all plots to have the same scaling, see the note below. Right now you will have to manually set the max value (sorry)
get_derivs_and_plot <- function(modobj,smooth_var,low_color=NULL,hi_color=NULL){
  this_font_size = font_size*1.25
  if (is.null(low_color)){low_color = "white"}
  if (is.null(hi_color)){hi_color = "firebrick"}
  derv<-derivatives(modobj,term=smooth_var, partial_match=T)
  derv<- derv %>%
    mutate(sig = !(0 >lower & 0 < upper))
  derv$sig_deriv = derv$derivative*derv$sig
  cat(sprintf("\nSig change: %1.2f - %1.2f\n",min(derv$data[derv$sig==T]),max(derv$data[derv$sig==T])))
  d1<- ggplot(data=derv) + geom_tile(aes(x = data, y = .5, fill = sig_deriv))
  
  # Set the gradient colors
  if (min(derv$derivative)>0) {
    d1 <- d1 + scale_fill_gradient(low = low_color, high = hi_color,limits = c(0,max(derv$derivative)))
    # If you want all plots to have the same scaling, this code can be used instead-- This is desirable if you have a factor-smooth model.
    #max_val = .01
    #d1 <- d1 +scale_fill_gradient(low = low_color,high = hi_color,limits = c(0,max_val),oob = squish)
  } else if (min(derv$derivative)<0 & max(derv$derivative)<0) {
    d1 <- d1 + scale_fill_gradient(low = hi_color, high = low_color,limits = c(min(derv$derivative),0))
  }else {
    d1 <- d1 +
      scale_fill_gradient2(low = "steelblue", midpoint = 0, mid = "white",
                           high = "firebrick",limits = c(min(derv$derivative),max(derv$derivative)))
    #this is only used to make the scales match on different plots
    #max_val = .01
    #d1 <- d1 +scale_fill_gradient(low = low_color,high = hi_color,limits = c(0,max_val),oob = squish)
  }
  
  d1 <- d1 + 
    labs(x = smooth_var,fill = sprintf("\u0394%s",smooth_var)) + 
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = this_font_size),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=this_font_size),
          legend.text = element_text(size = this_font_size),
          axis.title = element_text(size = this_font_size),
          legend.key.width = unit(1,"cm"),
          legend.position = "right",
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    guides(fill = guide_colorbar(reverse = F,direction = "horizontal",title.position = "top")) +
    geom_rect(aes(ymin=0,ymax=1,xmin=min(data),xmax=max(data)),color="black",fill="white",alpha = 0)
  return(d1)
}

## Parametric bootstrap of likelihood ratio test for nested models
## Takes your gamm model object as an input. If you are using `gam` this will have to be tweaked a little.
## Right now, this requires your term of interest to be the LAST term in the model.
pboot <- function(modelobj){
  library(lme4)
  numsims <- 10000 #This is the number of bootstrap simulations. This # should be higher in practice, probably something like 10,000
  
  df <- modelobj$gam$model #Get a data frame of all data used in the model (if using gam, this is modelobj$model)
  group_var_name <- tail(names(modelobj$lme$coefficients$random), n=1) # the name of the random effects variable
  thisResp <- as.character(modelobj$gam$terms[[2]])
  f1 <- modelobj$gam$formula #Get the formula (if using gam, this is modelobj$formula)
  theseVars <- attr(terms(f1),"term.labels") #Get the terms
  f2 <- reformulate(theseVars[0:(length(theseVars)-1)],response = thisResp) #Drop the last term from the formula. This is the simpler model now.
  
  # The bootstrap function takes an lme object, so we fit the models using gam to get the design matrix, and stick that matrix into lmer
  #Fit
  g1 <- gam(f1,data = df)
  g2 <- gam(f2,data = df)
  #print(summary(g1))
  #print(summary(g2))
  
  #Get the matrices
  mat1 <- model.matrix(g1)
  mat2 <- model.matrix(g2)
  
  #Tack on the response variable and grouping variable.
  group_var<- df[,group_var_name]
  #print(group_var)
  y <- df[,thisResp]
  #print(y)
  
  # Fit the models with `lmer`
  m1 <- lmer(y ~ -1 + mat1 + (1|group_var))
  m2 <- lmer(y ~ -1 + mat2 + (1|group_var))
  
  # Create the bootstrap distribution
  refdist <- PBrefdist(m1, m2, nsim=numsims) # note you can parallelize this (look at help for the function)
  pb <- PBmodcomp(m1, m2, ref = refdist) # Compare the models
  int_pval <- pb$test["PBtest","p.value"] # Pull out the p-value from the bootstrap test.
  
  # Now return the best model
  if (int_pval < .05) {
    pb$bestmod <- f1
  } else {
    pb$bestmod <- f2
  }
  return(pb)
}

#### FIT GAMM SMOOTH ####
##Function to fit a GAMM (measure ~ s(smooth_var, k = knots, fx = set_fx) + covariates, with a random intercept for participants with random=list(id_var =~0)) per each region in atlas and save out statistics
## Function adapted from Sydnor et al., 2023, https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/GAM_functions.R
gamm.fit.smooth <- function(measure, atlas, dataset, region, smooth_var, covariates, knots, id_var, random_slope=FALSE, set_fx = FALSE, stats_only = FALSE){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates))
  if (random_slope == FALSE){
    re_values=c(~1)
    re_names=id_var
    #gam.model <- gamm(modelformula, method = "REML",random = as.list(setNames(re_values,re_names)), data = gam.data)
  }
  if (random_slope == TRUE)
  {
    re_values=c(~1, as.formula(paste0("~0 +",smooth_var))) #the random effects have to be formulas
    re_names=c(id_var, id_var)
    #reStruct(as.list(setNames(values,names))) #this is how to create a random effects structure, see next line
  }
  gam.model <- gamm(modelformula, method = "REML",random = reStruct(as.list(setNames(re_values,re_names))), data = gam.data)
  gam.results <- summary(gam.model$gam)
  
  #GAMM derivatives
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  #Identify derivative significance window(s)
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv derivatives column where non-significant derivatives are set to 0
  
  #GAMM statistics
  #F value for the smooth term and GAM-based significance of the smooth term
  gam.smooth.F <- gam.results$s.table[3]
  gam.smooth.pvalue <- gam.results$s.table[4]
  
  #Calculate the magnitude and significance of the smooth term effect by comparing full and reduced models
  ##Compare a full model GAM (with the smooth term) to a nested, reduced model (with covariates only)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, covariates)) #no smooth term
  gam.nullmodel <- gamm(nullmodel, method = "REML",random = reStruct(as.list(setNames(re_values,re_names))), data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel$gam)
  
  ##Full versus reduced model anova p-value
  #We cannot use the ANOVA framework to compare GAMMs...because GAMMs are not truly wholly lme4 objects
  #Could substitute in the parametric bootstrap here if we wanted to, but it would take a while, and unlikely to be different
  #I uncommented this but it may need to be commented out
  #anova(gam.nullmodel$lme,gam.model$lme)
  #anova.smooth.pvalue <- anova.gam(gam.nullmodel$gam,gam.model$gam)$`p-value`[2]
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$gam$y - gam.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$gam$y - gam.nullmodel$gam$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  
  ### effect direction
  #I think this makes sense in the context of a regular continuous effect, but not in the case of an interaction in the next function
  mean.derivative <- mean(derv$derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    partialRsq <- partialRsq*-1}
  
  #Derivative-based temporal characteristics
  #Age of developmental change onset
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    change.onset <- min(derv$data[derv$sig==T])} #find first age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ #if gam derivative is never significant
    change.onset <- NA} #assign NA
  
  #Age of maximal developmental change
  if(sum(derv$sig) > 0){ 
    derv$abs_sig_deriv = round(abs(derv$sig_deriv),5) #absolute value significant derivatives
    maxval <- max(derv$abs_sig_deriv) #find the largest derivative
    window.peak.change <- derv$data[derv$abs_sig_deriv == maxval] #identify the age(s) at which the derivative is greatest in absolute magnitude
    peak.change <- mean(window.peak.change)} #identify the age of peak developmental change
  if(sum(derv$sig) == 0){ 
    peak.change <- NA}  
  
  #Age of decrease onset
  if(sum(derv$sig) > 0){ 
    decreasing.range <- derv$data[derv$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
    if(length(decreasing.range) > 0)
      decrease.onset <- min(decreasing.range) #find youngest age with a significant negative derivative
    if(length(decreasing.range) == 0)
      decrease.onset <- NA}
  if(sum(derv$sig) == 0){
    decrease.onset <- NA}  
  
  #Age of increase offset
  if(sum(derv$sig) > 0){ 
    increasing.range <- derv$data[derv$sig_deriv > 0] #identify all ages with a significant positive derivative (i.e., smooth_var indices where y is increasing)
    if(length(increasing.range) > 0)
      increase.offset <- max(increasing.range) #find oldest age with a significant positive derivative
    if(length(increasing.range) == 0)
      increase.offset <- NA}
  if(sum(derv$sig) == 0){ 
    increase.offset <- NA}  
  
  #Age of maturation
  if(sum(derv$sig) > 0){ 
    change.offset <- max(derv$data[derv$sig==T])} #find last age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ 
    change.offset <- NA}  
  
  full.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq, change.onset, peak.change, decrease.onset, increase.offset, change.offset)
  stats.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}

#### FIT GAMM SMOOTH x CONTINUOUS COVARIATE INTERACTION ####
##Function to fit a GAMM (measure ~ s(smooth_var, k = knots, fx = set_fx) + s(smooth_var, k = knots, fx = set_fx, by = int_var) + covariates)), with a random intercept for participants with random=list(id_var =~0)) per each region in atlas and save out statistics
## Function adapted from Sydnor et al., 2023, https://github.com/PennLINC/spatiotemp_dev_plasticity/blob/main/gam_models/GAM_functions.R
gamm.fit.smooth.int <- function(measure, atlas, dataset, region, smooth_var, int_var, id_var, covariates, knots, set_fx = FALSE, stats_only = TRUE, random_slope=FALSE){
  
  #Fit the gam
  dataname <- sprintf("%s.%s.%s", measure, atlas, dataset) 
  gam.data <- get(dataname)
  parcel <- region
  region <- str_replace(region, "-", ".")
  modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, k = %3$s, fx = %4$s) + s(%2$s, by = %5$s, k = %3$s, fx = %4$s) + %6$s", region, smooth_var, knots, set_fx, int_var, covariates))
  if (random_slope == FALSE){
    re_values=c(~1)
    re_names=id_var
    #gam.model <- gamm(modelformula, method = "REML",random = as.list(setNames(re_values,re_names)), data = gam.data)
  }
  if (random_slope == TRUE)
  {
    re_values=c(~1, as.formula(paste0("~0 +",smooth_var))) #the random effects have to be formulas
    re_names=c(id_var, id_var)
    #reStruct(as.list(setNames(values,names))) #this is how to create a random effects structure
  }
  gam.model <- gamm(modelformula, method = "REML",random = reStruct(as.list(setNames(re_values,re_names))), data = gam.data)
  gam.results <- summary(gam.model$gam)
  
  #GAMM derivatives
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(gam.model, term = sprintf('s(%s)',smooth_var), interval = "simultaneous", unconditional = F) #derivative at 200 indices of smooth_var with a simultaneous CI
  #Identify derivative significance window(s)
  derv <- derv %>% #add "sig" column (TRUE/FALSE) to derv
    mutate(sig = !(0 > lower & 0 < upper)) #derivative is sig if the lower CI is not < 0 while the upper CI is > 0 (i.e., when the CI does not include 0)
  derv$sig_deriv = derv$derivative*derv$sig #add "sig_deriv derivatives column where non-significant derivatives are set to 0
  
  #GAM statistics
  #F value for the smooth term and GAM-based significance of the smooth x interaction term
  gam.smooth.F <- gam.results$s.table[2,3]
  gam.smooth.pvalue <- gam.results$s.table[2,4]
  
  #Calculate the magnitude and significance of the smooth term effect by comparing full and reduced models
  ##Compare a full model GAMM (with the smooth term and the interaction term) to a nested, reduced model (with smooth term only)
  nullmodel <-   as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, covariates)) #no smooth term
  gam.nullmodel <- gamm(nullmodel, method = "REML",random = reStruct(as.list(setNames(re_values,re_names))), data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel$gam)
  
  ##Full versus reduced model anova p-value
  #We cannot use the ANOVA framework to compare GAMMs...because GAMMs are not truly wholly lme4 objects
  #Could substitute in the parametric bootstrap here if we wanted to, but it would take a while, and unlikely to be different
  #I uncommented this but it may need to be commented out
  # print(anova(gam.nullmodel$lme,gam.model$lme,test='Chisq'))
  # anova.smooth.pvalue <- anova.gam(gam.nullmodel$lme,gam.model$lme,test='Chisq')$`Pr(>Chi)`[2]
  
  ##Full versus reduced model direction-dependent partial R squared
  ### effect size
  sse.model <- sum((gam.model$gam$y - gam.model$gam$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$gam$y - gam.nullmodel$gam$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model)/sse.nullmodel
  #print(region);print(sse.model);print(sse.nullmodel);print(gam.smooth.pvalue)
  
  ### effect direction does not make sense in the context of an interaction, really, at least not for our purposes
  # mean.derivative <- mean(derv$derivative)
  # if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
  #   partialRsq <- partialRsq*-1}
  
  #Derivative-based temporal characteristics
  #Age of developmental change onset
  if(sum(derv$sig) > 0){ #if derivative is significant at at least 1 age
    change.onset <- min(derv$data[derv$sig==T])} #find first age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ #if gam derivative is never significant
    change.onset <- NA} #assign NA
  
  #Age of maximal developmental change
  if(sum(derv$sig) > 0){ 
    derv$abs_sig_deriv = round(abs(derv$sig_deriv),5) #absolute value significant derivatives
    maxval <- max(derv$abs_sig_deriv) #find the largest derivative
    window.peak.change <- derv$data[derv$abs_sig_deriv == maxval] #identify the age(s) at which the derivative is greatest in absolute magnitude
    peak.change <- mean(window.peak.change)} #identify the age of peak developmental change
  if(sum(derv$sig) == 0){ 
    peak.change <- NA}  
  
  #Age of decrease onset
  if(sum(derv$sig) > 0){ 
    decreasing.range <- derv$data[derv$sig_deriv < 0] #identify all ages with a significant negative derivative (i.e., smooth_var indices where y is decreasing)
    if(length(decreasing.range) > 0)
      decrease.onset <- min(decreasing.range) #find youngest age with a significant negative derivative
    if(length(decreasing.range) == 0)
      decrease.onset <- NA}
  if(sum(derv$sig) == 0){
    decrease.onset <- NA}  
  
  #Age of increase offset
  if(sum(derv$sig) > 0){ 
    increasing.range <- derv$data[derv$sig_deriv > 0] #identify all ages with a significant positive derivative (i.e., smooth_var indices where y is increasing)
    if(length(increasing.range) > 0)
      increase.offset <- max(increasing.range) #find oldest age with a significant positive derivative
    if(length(increasing.range) == 0)
      increase.offset <- NA}
  if(sum(derv$sig) == 0){ 
    increase.offset <- NA}  
  
  #Age of maturation
  if(sum(derv$sig) > 0){ 
    change.offset <- max(derv$data[derv$sig==T])} #find last age in the smooth where derivative is significant
  if(sum(derv$sig) == 0){ 
    change.offset <- NA}  
  
  full.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq, change.onset, peak.change, decrease.onset, increase.offset, change.offset)
  stats.results <- cbind(parcel, gam.smooth.F, gam.smooth.pvalue, partialRsq)
  if(stats_only == TRUE)
    return(stats.results)
  if(stats_only == FALSE)
    return(full.results)
}
