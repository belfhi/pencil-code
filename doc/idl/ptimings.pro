; $Id: ptimings.pro,v 1.4 2006-11-08 07:05:21 brandenb Exp $
;
red = 122
blue = 55
organge = 155
yellow = 188
fg = !p.color
if !d.name eq 'PS' then begin
  device,xsize=18,ysize=12,yoffset=3
  !p.charthick=3 & !p.thick=3 & !x.thick=3 & !y.thick=3
  ;red=0 & blue=0 & fg=0
end
;
;  mv idl.ps ../figs/ptimings.ps
;
fact=1.
a=rtable('timings.dat',2,head=1)
b=rtable('horseshoe.dat',2,head=1)
c=rtable('kabul.dat',2,head=1)
d=rtable('horseshoe_mega.dat',2,head=1)
e=rtable('pdc_06.dat',2,head=1)
e2=rtable('lenngren_06.dat',2,head=1)
f=rtable('steno_06.dat',2,head=1)
n1=reform(a(0,*)) & t1=reform(a(1,*))
n2=reform(b(0,*)) & t2=reform(b(1,*))
n3=reform(c(0,*)) & t3=reform(c(1,*))
n4=reform(d(0,*)) & t4=reform(d(1,*))
n5=reform(e(0,*)) & t5=reform(e(1,*))
n5b=reform(e2(0,*)) & t5b=reform(e2(1,*))
n6=reform(f(0,*)) & t6=reform(f(1,*))
;
save_state
sym = texsyms()
!p.multi=0
!p.charsize=1.4
!x.title='# of procs'
!y.title=sym.mu+'!6s/step/point'
!x.range=[.8,160]
!y.range=fact*[.05,15]
!x.style=3
!y.style=3

plot_oo, n1, fact*t1, PSYM=-1, LINE=fg
oplot,   n2, fact*t2, PSYM=-5,li=0, COLOR=red
oplot,   n3, fact*t3, PSYM=-6,li=2, COLOR=blue
oplot,   n4, fact*t4, PSYM=-7,li=3, COLOR=organge
oplot,   n5, fact*t5, PSYM=-4,li=0, COLOR=fg
oplot,   n5b, fact*t5b, PSYM=-4,li=0, COLOR=fg
oplot,   n6, fact*t6, PSYM=-5,li=2, COLOR=red
;
esrg_legend, ['!6Origin3000!X', '!6Horseshoe!X', '!6KIS cluster!X', 'GigaBits!X', '!6Steno', '!6Luci/Lenn'], $
    LINE=[1,0,2,3,2,0], $
    COLOR=[fg,red,blue,organge,red,fg], $
    PSYM=[-1,-5,-6,-7,-5,-4], $
    SPOS='tr', /BOX

restore_state
;
; xx=[1,120] & oplot,xx,4./xx^.7
print,'import ptimings.jpg'
print,'scp2 ptimings.jpg $scr/ccp2001'
;
end
