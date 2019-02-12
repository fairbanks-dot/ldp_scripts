import os

os.chdir('K:/18_Projects/Shishmaref AP Erosion Control/2-Survey/1-DataRaw/Photogrammetry/SHISHMAREF/SHISHMAREF_ISLAND_3D/101MEDIA')

# K:\18_Projects\Shishmaref AP Erosion Control\2-Survey\1-DataRaw\Photogrammetry\SHISHMAREF\SHISHMAREF_ISLAND_3D\101MEDIA
# this will check if you are in the directory you think you are
#print(os.getcwd())

# c is the first number you want to ADD for rename the file
c = 2000

#print('Re-naming files in: ', os.getcwd())
for f in os.listdir():
    #print(f)
    # this if statement will check if file found is NOT Thumbs.db and if has '.' in it, so it wont be a folder
    if f != 'Thumbs.db' and f.find('.') > 0:
        f_name, f_ext = os.path.splitext(f)
        #prefix, suffix = f_name.split('_')
        #n = int(suffix) + c
        #new_name = '{}_{}{}'.format(prefix, n,f_ext)
        #os.rename(f, new_name)
        
        new_name = '{}{}'.format(c, f_ext)
        
        #print(new_name)
        os.rename(f, new_name)
        c = c + 1
        #print(new_name)