pro idl_1

;#### Set grid ####

  numr=130
  numz=130
  numphi=256
  pi=3.1415
  omega=0.112


;#### Create arrays of input file names #### 

  file = 'fileslist'
  OPENR, lun, file, /GET_LUN
  filename = ''
  line = ''
  WHILE NOT EOF(lun) DO BEGIN & $
    READF, lun, line & $
    filename = [filename, line] & $
  ENDWHILE
  FREE_LUN, lun

  numfile=size(filename,/n_elements)-1
;  print,size(filename,/n_elements)
;  print,filename


;#### Set up arrays for plots ####

  timearr=dblarr(numfile)
  rhomaxarr=dblarr(numfile)
  rho=dblarr(numr,numz,numphi)
  rho_eq=dblarr(numphi,numr)
  R=dblarr(numr)
  theta=dblarr(numphi)

  for i=0,numfile-1 do begin
    timearr[i]=2.0*pi/omega * i/20.0 
  endfor  
  
  for i=0,numr-1 do begin
    R[i]=i
  endfor

  for i=0,numphi-1 do begin
    theta[i]=i*2*pi/(numphi-1.0) 
  endfor

;for i=0,numfile-1 do begin
;  print, filename[i]
;endfor


;  print, timearr

;#### Plotting commands ####

thick = 1
!p.thick = thick
!p.charthick = thick
!x.thick = thick
!y.thick = thick


device, decomposed=0
loadct,5


;#### Loop over files ####

for k=1,numfile do begin
  print, filename[k]

;#### Read density ####

  openr,1,filename[k],/f77_unformatted,/swap_endian
  readu,1,rho
  close,1

;#### Print vertical slice ####  

  cc=indgen(210)+10

  
;#### Print horizontal slice ####  
  for i=0,numphi-1 do begin
    for j=0,numr-1 do begin
      rho_eq[i,j]=rho[j,2,i]   
    endfor
  endfor

  
  print, "plotnow"
  polar_contour, rho_eq^0.25, theta, R, /DITHER, /fill, nlevels=210,$
  c_colors=cc, xrange=[-numr,numr],yrange=[-numr,numr],/xstyle,/ystyle  
  
;  cc=indgen(210)+40
  

;#### Max rho ####  
  rhomaxarr[k-1]=max(rho)
  
endfor


;#### Print max density as function of time ####
  device,filename='den.eps',/encapsulated,/portrait,bit=24,/color 
  DEVICE, XSIZE=5.0, YSIZE=4.5, /INCHES
  plot,timearr,rhomaxarr, yrange=[0,1.1],/ystyle
  device,/close
  
;  print,rhomaxarr

  set_plot,'x'
  
;#### End of code ####  
  
  
;  p=plot(timearr,rhomaxarr)

;  print,p
  
;for i=0,numphi-1 do begin
;  for j=0,numr-1 do begin
;    rho_eq[i,j]=rho[j,73,i]   
;  endfor
;endfor

;cc=indgen(210)+40


;polar_contour, rho_eq^0.2, theta, R, /DITHER, /fill, nlevels=210,$
;	c_colors=cc

	
	
end
