git remote add origin https://github.com/swansonalan/gwRaster
git remote set-url origin git@github.com:swansonalan/gwRaster.git

cd /home/aswanson/code/gwRaster
gid add .
git commit -m "initial commit"
git push -u origin master
git pull



# ssh code~~~~~~~~~~~~
# ssh-keygen -t ed25519 -C "mtskimtb@gmail.com"
# eval "$(ssh-agent -s)"
# ssh-add ~/.ssh/id_ed25519c
#

# package development code ~~~

usethis::use_build_ignore("ignore")
usethis::use_package("geodist")
usethis::use_package("raster")
usethis::use_package("abind")
usethis::use_package("foreach")
usethis::use_package("doParallel")
usethis::use_mit_license("Alan Swanson")



library(roxygen2); # Read in the roxygen2 R package
roxygenise();      # Builds the help files
devtools::document()
devtools::install()

