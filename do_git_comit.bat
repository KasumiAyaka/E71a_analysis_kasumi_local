if "%1"=="" goto :usage


git add .

git commit -m "%1"

git push



:usage

enter git commit message!
