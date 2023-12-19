;test initial data
wdir='h:\red_data.pol\Mkn1040_161124\'
ima=readfits(wdir+'s14720206.fts')
y_c=center_frames_WOLL2(ima)
print,y_c
end