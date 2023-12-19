objdir='Mrk1018'
dir='/data6/SCORPIO/sppol_pipeline_v2023.8/'+objdir+'/'
h=headfits(dir+'avg_spectra.fit')
Nx=sxpar(h,'NAXIS1')
lam=indgen(Nx)*sxpar(h,'CDELT1')+sxpar(h,'CRVAL1')
cdel=sxpar(h,'CDELT1') & crva=sxpar(h,'CRVAL1')
lam_st=4800 & lam_fin=7800											;Wavelength range
z=0.0255
;emlines=[6563,4861,4341]											;Broad lines to oplot //[6563,4861,4341,1549,2798,1909]
;narlines=[6550,6585,5007,4959]										;Narrow lines to oplot //[6550,6585,5007,4959]
plot_atm = 0														;Mark atm. abs. bands: 0 - NO, 1 - YES
;atm = [[6860,6917]]										;Atmosperic bands *OPTIONAL* H20=[7150,7350] / B=[6860,6917] / A=[7590,7720]
numexp=[sxpar(h,'NUMEXP1'),sxpar(h,'NUMEXP2'),sxpar(h,'NUMEXP3')]
date_obs=sxpar(h,'DATE-OBS')
ww=[2,2,2] 															;Binning windows for OBJ, NON-POL std, POL std
pa=[sxpar(h,'PA1'),sxpar(h,'PA2'),sxpar(h,'PA3')]
names=[sxpar(h,'NAME1'),sxpar(h,'NAME2'),sxpar(h,'NAME3')]
chr=2.6


;===============================================================

R=where(lam gt lam_st and lam lt lam_fin,ind)

;Channel transmission Dq & Du from NON-POL
avg_spectra=readfits(dir+'avg_spectra.fit')
cgdisplay, wid=0
!p.multi=[0,1,2]
stnum=1
cgplot, lam(R),avg_spectra(R,0,stnum)/avg_spectra(R,1,stnum), yrange=[median(avg_spectra(R,0,stnum)/avg_spectra(R,1,stnum))-0.2,median(avg_spectra(R,0,stnum)/avg_spectra(R,1,stnum))+0.2], $
		xrange=[lam_st-10,lam_fin+10]
		Dq=lowess(indgen(ind),avg_spectra(R,0,stnum)/avg_spectra(R,1,stnum),ind/2,2,1)
		cgoplot,lam(R),Dq,color='gold',thick=3
cgplot, lam(R),avg_spectra(R,2,stnum)/avg_spectra(R,3,stnum), yrange=[median(avg_spectra(R,2,stnum)/avg_spectra(R,3,stnum))-0.2,median(avg_spectra(R,2,stnum)/avg_spectra(R,3,stnum))+0.2], $
		 xrange=[lam_st-10,lam_fin+10]
		Du=lowess(indgen(ind),avg_spectra(R,2,stnum)/avg_spectra(R,3,stnum),ind/6,1,1)
		cgoplot,lam(R),Du,color='gold',thick=3
writefits, dir+'Dq.fts',dq,h
print,'Dq.fts is stored! ', dir
writefits, dir+'Du.fts',du,h
print,'Du.fts is stored! ', dir

;Stokes plots
tot_spec=fltarr(ind,3)
q=fltarr(ind,max(numexp),3) & u=fltarr(ind,max(numexp),3)
spectra=readfits(dir+'spectra.fit')

start_ps, dir+'stokes_'+objdir+'.eps'
for ob=0,2 do begin
	print, '----', names(ob), '----'
	;Calculate I
	tot_spec(*,ob)=avg_spectra(R,0,ob)+avg_spectra(R,1,ob)*Dq+avg_spectra(R,2,ob)+avg_spectra(R,3,ob)*Du
	;Correct depolarization in the cubes
	avg_ratio=fltarr(ind,2,numexp(ob))
	for ex=0,numexp(ob)-1 do begin
		for j=0,1 do begin
			for kx=0,ind-1 do begin
				avg_ratio(kx,j,ex)=spectra(R(kx),j*2,ex,ob)/spectra(R(kx),j*2+1,ex,ob)
			endfor
		endfor
	endfor
	for kx=0,ind-1 do begin
		for j=0,1 do begin
			robomean,avg_ratio(kx,j,*),1,0.5,mean
			avg_ratio(kx,j,*)=avg_ratio(kx,j,*)/median(avg_ratio(kx,j,*))
		endfor
	endfor
	avgd=fltarr(2,numexp(ob))
		for j=0,1 do begin
			for ke=0,numexp(ob)-1 do begin
				avgd(j,ke)=median(avg_ratio(*,j,ke))
			endfor
			robomean,avgd(j,*),3,0.5,mean
			avgd(j,*)=avgd(j,*)/median(avgd(j,*))
		endfor
	;Calculates Stokes Q and U
	nexp=numexp(ob)
	for ex=0,nexp-1 do begin
		q(*,ex,ob)=-(spectra(R,2,ex,ob)-spectra(R,3,ex,ob)*Du*avgd(1,ex)) / (spectra(R,2,ex,ob)+spectra(R,3,ex,ob)*Du*avgd(1,ex))
		u(*,ex,ob)=(spectra(R,0,ex,ob)-spectra(R,1,ex,ob)*Dq*avgd(0,ex)) / (spectra(R,0,ex,ob)+spectra(R,1,ex,ob)*Dq*avgd(0,ex))
	endfor
	;Rotate to celestial plane
		qn=q*cos(-2*pa(ob)*!pi/180.0)+u*sin(-2*pa(ob)*!pi/180.0)
		un=u*cos(-2*pa(ob)*!pi/180.0)-q*sin(-2*pa(ob)*!pi/180.0)
		q=qn & u=un
	;Stokes binning
	;:построение узлов
	w = ww(ob) & Nx=ind
	xpos=fltarr(Nx)
		Npos=Nx/w-1
			if ob eq 0 then	polarization=fltarr(Npos,6,3)
		xpos(0:Npos-1)=findgen(Npos)*w+w/2
	Qres=fltarr(Npos) & Ures=fltarr(Npos)
	d_Qres=fltarr(Npos) & d_Ures=fltarr(Npos)
		;:робастная оценка среднего и rms в окне W рх
		for k=0,Npos-1 do begin
			if w ne 1 then begin
				robomean,Q(xpos(k)-w/2:xpos(k)+w/2-1,0:numexp(ob)-1,ob),1,0.5,mean,rms
				Qres(k)=mean & d_Qres(k)=rms;/3
				robomean,U(xpos(k)-w/2:xpos(k)+w/2-1,0:numexp(ob)-1,ob),1,0.5,mean,rms
				Ures(k)=mean & d_Ures(k)=rms;/3
			endif else begin
				robomean,Q(k,0:numexp(ob)-1,ob),1,0.5,mean,rms
				Qres(k)=mean & d_Qres(k)=rms;/3
				robomean,U(k,0:numexp(ob)-1,ob),1,0.5,mean,rms
				Ures(k)=mean & d_Ures(k)=rms;/3
			endelse
		endfor

tot_spec(*,ob)=tot_spec(*,ob)/1.0e4
cgdisplay, wid=1, xsize=1350, ysize=1500
	!p.multi=[0,1,5]
	!p.charsize=1.8

;I panel
	cgplot, lam(R), tot_spec(*,ob),xrange=[lam_st,lam_fin], ytitle='F, 10$\up4$ ADU', xtickformat="(A1)", charsize=chr, $
		title=names(ob)+'   '+date_obs+'   BTA+SCORPIO-2', font=2, yrange=[0,max(tot_spec(*,ob))*1.1], pos=[0.11,0.09+0.172*4,0.95,0.09+0.172*5], $
		yTICKINTERVAL=10.0
			for ll=0,N_elements(emlines)-1 do cgoplot, [1,1]*emlines(ll)*(1+z), [-1e5,1.2*max(tot_spec(*,ob))], color='blue',linestyle=3
			for ll=0,N_elements(narlines)-1 do cgoplot, [1,1]*narlines(ll)*(1+z), [-1e5,1.2*max(tot_spec(*,ob))], linestyle=2
				if plot_atm eq 1 then begin
					for ia=0,N_elements(atm(0,*))-1 do begin
						vel=atm(*,ia)
						cgPolygon, [vel(0), vel(1),  vel(1),  vel(0), vel(0)], [0,0, max(tot_spec(*,ob))*1.09, max(tot_spec(*,ob))*1.09, 0], $
           					 fCOLOR='blk3', /overplot, /fill, color='blk3'
					endfor
				endif
        cgoplot, lam(R), tot_spec(*,ob)
			cgtext,lam_fin-0.052*(lam_fin-lam_st), max(tot_spec(*,ob))*0.9, '(1)', font=2, charsize=1.1

;Q panel
	cgplot, xpos(0:Npos-1)*cdel+lam_st, qres*100, yrange=[median(qres*100)-14.5,median(qres*100)+14.5], xrange=[lam_st,lam_fin], ytitle='Q, %', psym=10, $
		err_ylow=d_Qres*100, err_yhigh=d_Qres*100, err_width=0, font=2, pos=[0.11,0.09+0.172*3,0.95,0.09+0.172*4], xtickformat="(A1)", charsize=chr, yminor=4
			for ll=0,N_elements(emlines)-1 do cgoplot, [1,1]*emlines(ll)*(1+z), [-10,10], color='blue',linestyle=3
			for ll=0,N_elements(narlines)-1 do cgoplot, [1,1]*narlines(ll)*(1+z), [-10,10], linestyle=2
				if plot_atm eq 1 then begin
					mi=median(qres*100)-2.45 & ma=median(qres*100)+2.45
					for ia=0,N_elements(atm(0,*))-1 do begin
						vel=atm(*,ia)
						cgPolygon, [vel(0), vel(1),  vel(1),  vel(0), vel(0)], [mi,mi,ma,ma,mi], $
           					 fCOLOR='blk3', /overplot, /fill, color='blk3'
					endfor
				endif
        cgoplot, xpos(0:Npos-1)*cdel+lam_st, qres*100, psym=10,err_ylow=d_Qres*100, err_yhigh=d_Qres*100, err_width=0
      	  cgtext,lam_fin-0.052*(lam_fin-lam_st), (median(qres*100)+1.5), '(2)', font=2, charsize=1.1
				robomean,Qres*100,3,0.5,av,rms
				print, 'Mean Q = ', av,'+-', rms, '%'

;U panel
	cgplot, xpos(0:Npos-1)*cdel+lam_st, ures*100,  yrange=[median(ures*100)-14.5,median(ures*100)+14.5]+0.2, xrange=[lam_st,lam_fin], ytitle='U, %', psym=10, $
		err_ylow=d_Ures*100, err_yhigh=d_Ures*100, err_width=0, font=2, pos=[0.11,0.09+0.172*2,0.95,0.09+0.172*3], xtickformat="(A1)", charsize=chr, yminor=4
			for ll=0,N_elements(emlines)-1 do cgoplot, [1,1]*emlines(ll)*(1+z), [-10,10], color='blue',linestyle=3
			for ll=0,N_elements(narlines)-1 do cgoplot, [1,1]*narlines(ll)*(1+z), [-10,10], linestyle=2
				if plot_atm eq 1 then begin
					mi=median(ures*100)-2.45+0.2 & ma=median(ures*100)+2.45+0.2
					for ia=0,N_elements(atm(0,*))-1 do begin
						vel=atm(*,ia)
						cgPolygon, [vel(0), vel(1),  vel(1),  vel(0), vel(0)], [mi,mi,ma,ma,mi], $
           					 fCOLOR='blk3', /overplot, /fill, color='blk3'
					endfor
				endif
        cgoplot, xpos(0:Npos-1)*cdel+lam_st, ures*100, psym=10,err_ylow=d_ures*100, err_yhigh=d_ures*100, err_width=0
      	  cgtext,lam_fin-0.052*(lam_fin-lam_st), (median(ures*100)+1.5)+0.2, '(3)', font=2, charsize=1.1
				robomean,Ures*100,3,0.5,av,rms
				print, 'Mean U = ', av,'+-', rms, '%'

;PD panel
	p=sqrt(qres*qres + ures*ures)
	dp=sqrt((qres*d_qres)^2 + (ures*d_ures)^2)/p
	cgplot, xpos(0:Npos-1)*cdel+lam_st, p*100,  yrange=[median(P*100)-14.5,median(P*100)+14.5], xrange=[lam_st,lam_fin], ytitle='P, %', psym=10, $
		err_ylow=dp*100, err_yhigh=dp*100, err_width=0, font=2, pos=[0.11,0.09+0.172*1,0.95,0.09+0.172*2], xtickformat="(A1)", charsize=chr, yminor=4
			for ll=0,N_elements(emlines)-1 do cgoplot, [1,1]*emlines(ll)*(1+z), [-10,10], color='blue',linestyle=3
			for ll=0,N_elements(narlines)-1 do cgoplot, [1,1]*narlines(ll)*(1+z), [-10,10], linestyle=2
				if plot_atm eq 1 then begin
					mi=median(p*100)-2.45 & ma=median(p*100)+2.45
					for ia=0,N_elements(atm(0,*))-1 do begin
						vel=atm(*,ia)
						cgPolygon, [vel(0), vel(1),  vel(1),  vel(0), vel(0)], [mi,mi,ma,ma,mi], $
           					 fCOLOR='blk3', /overplot, /fill, color='blk3'
					endfor
				endif
        cgoplot, xpos(0:Npos-1)*cdel+lam_st, p*100, psym=10,err_ylow=dp*100, err_yhigh=dp*100, err_width=0
       		cgtext,lam_fin-0.052*(lam_fin-lam_st), (median(p*100)+1.5), '(4)', font=2, charsize=1.1
				robomean,P*100,3,0.5,av,rms
				print, 'Mean PD = ', av,'+-', rms, '%'

;PA panel
	phi=p
	for i=0,npos-1 do begin
		qtmp=qres(i) & 	utmp=ures(i)
        phi(i)=atan(qtmp,utmp)
	;	phi(i)=arctan(qtmp,utmp)
	endfor
	;Correct 180 deg if needed
	;	RF=WHERE(phi gt 150, indp) & if indp gt 0 then phi(RF)=phi(RF)-180
	;	RF=WHERE(phi lt 60, indp) & if indp gt 0 then phi(RF)=phi(RF)+180
	letter="152B
	dphi=28.65*dp/p
	cgplot, xpos(0:Npos-1)*cdel+lam_st, phi,  xrange=[lam_st,lam_fin], charsize=chr, ytitle='!9' + String(letter) + '!X'+', deg', psym=10, $
		err_ylow=dphi, err_yhigh=dphi, err_width=0, yrange=[median(phi)-80,median(phi)+80],/ERR_CLIP, font=2, $
		xtitle='Wavelength, '+cgsymbol('Angstrom'), pos=[0.11,0.09+0.172*0,0.95,0.09+0.172*1]
			for ll=0,N_elements(emlines)-1 do cgoplot, [1,1]*emlines(ll)*(1+z), [median(phi)-85,median(phi)+85], color='blue',linestyle=3
			for ll=0,N_elements(narlines)-1 do cgoplot, [1,1]*narlines(ll)*(1+z), [median(phi)-85,median(phi)+85], linestyle=2
				if plot_atm eq 1 then begin
					mi=median(phi)-78 & ma=median(phi)+78
					for ia=0,N_elements(atm(0,*))-1 do begin
						vel=atm(*,ia)
						cgPolygon, [vel(0), vel(1),  vel(1),  vel(0), vel(0)], [mi,mi,ma,ma,mi], $
           					 fCOLOR='blk3', /overplot, /fill, color='blk3'
					endfor
				endif
        cgoplot, xpos(0:Npos-1)*cdel+lam_st, phi, psym=10,err_ylow=dphi, err_yhigh=dphi, err_width=0, /err_clip
       		cgtext,lam_fin-0.052*(lam_fin-lam_st), (median(phi)+55), '(5)', font=2, charsize=1.1
				robomean,phi,3,0.5,av,rms
				print, 'Mean PHI = ', av,'+-', rms, '%'

	if ob ne 2 then  ERASE

;Write OBJ to txt
file=dir+objdir+'.txt'
	if ob eq 0 then begin
		openw, 4, file, width=200
		for k=0,Npos-1 do begin
			printf, 4, xpos(k)*cdel+lam_st, qres(k)*100, d_qres(k)*100, ures(k)*100, d_ures(k)*100, p(k)*100, dp(k)*100, phi(k), dphi(k)
		endfor
	close, 4
	endif

endfor

stop_ps

print, '*****Plot in eps was written to ',dir+'stokes_'+objdir+'_TEST.eps', ' (3 pages)'
print, '*****OBJ specpol was written to '+dir+objdir+'.txt', ' (wavelength,Q,dQ,U,dU,P,dP,PHI,dPHI)'

end
