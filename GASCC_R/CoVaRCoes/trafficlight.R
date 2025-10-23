trafficlight <- function(data, name, plottitle) {
  colnames(data) <- c("A&F","C","P","NFM","Che","S","E","LM","U","EP","T","BM","EE","ME","HA","F&D","T&A","Bio",
                      "RE","CT","L&S","Com","B&D","ND","Comp","Med","Comm","B","NBF","A","Beau","or")
  rownames(data) <- name
  if (max(data[,1:31])==3 ){
    col<-c( "#41Ab5d", "#FFFF33","#FF8C00")
    
    p <- superheat(data[,1:31], 
                  left.label.size = 0.33,
                  bottom.label.size = 0.25, 
                  scale=F,
                  order.rows = order(data$or),
                  heat.pal=col,
                  left.label.text.col="#2F4F4F", 
                  left.label.text.angle = 0,    
                  left.label.text.alignment = "left",    
                  left.label.col=c("#F7FBFF","#DEEBF7"),
                  bottom.label.col=c("#F0FFFF"),  
                  bottom.label.text.col="#2F4F4F",   
                  bottom.label.text.angle=80,    
                  bottom.label.text.alignment = "center", 
                  bottom.label.text.size = 3.5, left.label.text.size = 3.8,
                  title = (plottitle),    
                  title.size=7,  
                  grid.hline.col = "#525252",      
                  grid.hline.size = 0.5,        
                  grid.vline.col = "#525252",     
                  grid.vline.size =0.5,
                  legend=F,
                  force.grid.vline=T,force.grid.hline=T,
                  left.label = 'variable',
                  heat.pal.values = c(0, 0.5, 1)
    )
  }else{
    col<-c( "#41Ab5d", "#FFFF33")
    p <- superheat(data[,1:31], 
                  left.label.size = 0.33,
                  bottom.label.size = 0.25, 
                  scale=F,
                  order.rows = order(data$or),
                  heat.pal=col,
                  left.label.text.col="#2F4F4F", 
                  left.label.text.angle = 0,    
                  left.label.text.alignment = "left",    
                  left.label.col=c("#F7FBFF","#DEEBF7"),
                  bottom.label.col=c("#F0FFFF"),  
                  bottom.label.text.col="#2F4F4F",   
                  bottom.label.text.angle=80,    
                  bottom.label.text.alignment = "center", 
                  bottom.label.text.size = 3.5, left.label.text.size = 3.8,
                  title = (plottitle),    
                  title.size=7,  
                  grid.hline.col = "#525252",      
                  grid.hline.size = 0.5,        
                  grid.vline.col = "#525252",     
                  grid.vline.size =0.5,
                  legend=F,
                  force.grid.vline=T,force.grid.hline=T,
                  left.label = 'variable',
                  heat.pal.values = c(0, 0.5, 1)
    )
  }
  
  return(p)
}