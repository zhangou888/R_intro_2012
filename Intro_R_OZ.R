## Introduction to R and R Programming
## Ou Zhang 09/26/2012
## Please create a folder in your C-drive called 'Rtest'
## you can also operate this function by 
dir.create('C:/Rtest') 

# please save all the files in the 'Rtest' folder.
# Please be aware that "R is case sensitive"!

## two types of help command:
# option 1:
help(vector)
?vector

# option 2:
help.search("vector")
??Matrix

# if you know an R function but you don't know keywords inside this function, 
# you just press 'tab' key in RStudio for help.
read.csv()  # please press 'tab'.

### Quit R
q()
quit()

## demonstration demo().
demo()
demo(image)
demo(graphics)

## get work directory.
getwd()

## set work directory.
setwd("c:/Rtest")

# To list all of the objects that are in memory, use 'ls()'
ls()
# To remove objects, use 'rm'
rm()

#----------------------------- R basic ------------------------------
#--- Kill current operation: Esc or Ctrl-C


#----------------------------- R as a calculator -----------
# --- basic operators
1+3  # addition
2-1  # subtraction
5*6  # multiplication
20/4 # division
5^2  # power

# --- built in functions
# general function call in R: function_call(arguments)
# e.g.

exp(5)               # exponential function
log(5)               # logarithm
cos(0)               # cosine
sum(c(1,2,3,4,5,6))  # sum
mean(c(1,2,3,4,5,6)) # mean
sd(c(1,2,3,4,5,6))   # standard deviation

# functions can be nested; 
# "from the inside out" according to the following order of operations
# 1. terms inside parentheses or brackets
# 2. exponents and roots
# 3. multiplication and division As they appear left to right
# 4. addition and subtraction As they appear left to right

exp(mean(c(1,2,3,4,5,6)))^5
# 1. calcuate mean of this number bundle.
# 2. take the exp() of the mean value
# 3. take that to the 5th power

# more examples below.
4+6*2
(4+6)*2
4*6+2
(4*6)+2
4*(6+2)
4*(6+2)^2
(4*(6+2))^2

## ---------- Objects and data types in R --------------------------------
#--- Assignments to variables --- 
x <- 5                         # store 5 in variable (object) x "x gets 5"
x                              # print the value of x on screen
(x <- 5)                       # assign and print -> use parentheses
(y <- 3)                       # store 3 in variable y
(z <- x + y)                   # sum of two variables option1
(z <- sum(x, y))               # sum of two variables option2
(z <- x*y)
(z <- x^y)
(z <- sqrt(x))
(z <- exp(y))
(z <- log(exp(y)))

#--- Vectors ---                            
# ordered set of elements of the same type
# create vector from 1 to 5
(vec0 <- 1:5)    
# create a vector using c() combine
(vec1 <- c(1, 2, 3))
#display the structure of an R object VERY USEFUL!
str(vec1)
# repeat elements
(vec2 <- rep(0, 5)) 
# generate sequence
(vec3 <- seq(from = 0, to = 1, by = 0.2))
#sample random numbers with replacement
(vec4 <- sample(1:10, 20, replace = TRUE)) 
# Vector of characters, characters are set between inverted commas " " or ' '
(vec5 <- c("a", "b", "c"))
(vec6 <- c(3, "a", 2, "b"))                        
str(vec5) 
str(vec6)
#combine vectors
(vec7 <- c(vec1,vec0)) 

# some vector operations that are useful
#length of a vector
length(vec0)   
# reverse vector order
rev(vec0)        
# logical function of vector
is.vector(vec0)
is.vector(NULL)

#---- Matrices ---                        
#2-dimensional array of elements of the same type
mat1 <- matrix(1:10, ncol = 2)                 #elements column-wise
mat1 <- matrix(1:10, ncol = 2, byrow = TRUE)   #elements row-wise
mat2 <- t(mat1)                                #transpose

# matrix multiplication and hadamard product. check the difference
mat1%*%mat2                                    #matrix multiplication
mat1*mat1                                      #Hadamard product

#select elements                               
mat1[1,2]                                      #select element (1,2)
mat1[1,]                                       #select row 1 and all columns
mat1[,2]                                       #select column 2 and all rows
mat1[,-2]                                      #select all columns but column 2

#combine vectors to matrix
vec8 <- 1:5
vec9 <- 6:10

(mat3 <- rbind(vec8,vec9)) #combine vectors to a matrix, with vectors in rows
(mat3 <- cbind(vec8,vec9)) #combine vectors to a matrix, with vectors in columns
mat3                       #cbind/rbind recycle row and column names

#useful matrix functions
dim(mat1)                  #dimensions of array
dim(mat1)[1]               #number of rows
dim(mat1)[2]               #number of columns
ncol(mat1)                 #number of rows  
nrow(mat1)                 #number of columns
length(mat1)               #number of ELEMENTS in matrix (rows x columns)
is.matrix(mat1)

rowSums(mat1)              #sum of row
colSums(mat1)              #sum of column
rowMeans(mat1)             #mean of row
colMeans(mat1)             #mean of column

#naming matrices
rownames(mat1)                              #check rownames
(rownames(mat1) <- paste("row",1:5))        #set rownames
colnames(mat1) <- c("col1","col2")          #set colnames
rownames(mat1)
mat1

#--- Factors
#factors are a data type for categorical variables
# R will then treat it differently

(vec12 <- c(1,2,3,3,3,2,2,1,3,2,1))       #numeric vector
is.numeric(vec12)                         
is.factor(vec12)                          #it is not a factor  

f.vec12 <- factor(vec12)                  #specify vec12 as factor
is.factor(f.vec12)                        #it is a factor 

#specify vec12 as ordered factor (ordinal)
of.vec12 <- factor(vec12,ordered=TRUE)    

#specify different labels when coercing 
#(from 1,2,3 to a,b,c) 
f.vec12 <- factor(vec12,labels=c("a","b","c"))

#--- Data Frame ---                            
# similar to matrix, but different types of variables are possible 
# e.g. characters and numbers and boolean
# standard data format if you read data into R
df1 <- letters[1:10]                    # letters a...j
df2 <- 1:10                             # numbers 1 to 10
(df3 <- rep(c(TRUE,FALSE),5))           # a vector of TRUE and FALSE alternating
str(df1)                                       
str(df2)
str(df3)

#concatenate as data frame
(df11 <- data.frame(df2, df3))                   

#a data frame is not a matrix by this test!
is.matrix(df11)        
#a data frame is a data frame however
is.data.frame(df11)                            
df12 <- data.frame(df1, df2)                   
str(df12)
df12[5,]                                       #select 5th row

#--- Lists ---    
# ordered set of components that may have different (data) types ; 
# most flexible data structure
list1 <- list(df1, df2, vec1, mat1)     #list of letters, vectors, matrix
list1[[1]]                              #select first list element (letters)
list1[[1]][5]                           #pick out 5th letter

list1[[4]]                              #select matrix
list1[[4]][4,2]                         #pick out element (4,2)

names(list1) <- c("let","seq","vec","mat")  #label list
list1
list1$mat                                   #can be used instead of list1[[4]] 
list1$mat[,2]                               #select 2nd column

#--- Missing values - 'NA' 
x <- 1:10
x[c(1,2,5)] <- NA                    #Insert Missing Values on the places 1, 2 and 5
x                                    #x with NAs
na.omit(x)                           #x is returned without NAs; additionally: Position of the missings
na.pass(x)                           #x is returned without changes
na.fail(x)                           #error if x includes NAs (x is not returned)

sum(x)                               #Statistical functions in R may require that all values are nonmissing!
mean(x)                           
sum(x, na.rm = TRUE)                 #Remove missing values before the computation

is.na(x)                             #logical condition for missing values in x
any(is.na(x) == TRUE)                #Any missing values in x?
which(is.na(x) == TRUE)              #Position of the missing values in x

1/0                                  #returned value is infinite (Inf)
0/0                                  #returned value in not a number (NaN)
x[2] <- 0/0
is.na(x)                             #no distinction between NA and NaN/Inf
is.nan(x)                            #to distinguish 'ordinary' NAs from missing values that result from
#certain computation

#--------------------------------- end objects & data types -----------------------
#-------------------------- getting data into R ---------------------------------
#--- Save and Load
save(data,file="data.rda")                  #saves data file in a compressed format; use extensions like .rda, .Rda, .RData etc.
load("data.rda")                            #loads a data set in R format

#--- import/export foreign data files
# A package is needed - 'foreign' 
# please install these this package

# load library
library(foreign)                                     #import/export package

# read txt-file
data1 <- read.table("C:/Rtest/data1.txt", header = TRUE)  
head(data1)   # check first 6 rows
tail(data1)   # check last 6 rows
getwd()       # In R it is needed to change the directory! Otherwise quote path.
?read.table

# read spss files with and without labels
Dats1 <- read.spss(file = "hs0.sav", to.data.frame = TRUE)   
Dats2 <- read.spss(file = "hs1.sav", to.data.frame = TRUE, use.value.labels = FALSE)

head(Dats1)                                  #Column names and first 6 rows of the dataset
head(Dats2)
tail(Dats1)                                  #Column names and last 6 rows of the dataset

# read excel file .xls file in the workspace
# Theoretically, we can read .xls files in the workspace through these 
# two packages-'xlsx' or 'xlsReadWrite"
# But, in practice, it's always painful to do that. 
# A lot of literature, recommend convert .xls file to .csv file
# I am going to skip this section.
# loading package.
library(xlsReadWrite)
library(xlsx)
library(gdata)

# if you have .csv files, you can use read.csv() or read.table() to import data
# read.table()
data2 <- read.table("hs0.csv",header=T,sep=",")
# read.csv()
data3 <- read.csv("hs0.csv",header=T,sep=",")


# read SAS data file.
# need package-"foreign"
tb1 <- read.ssd("C:/Rtest", "hs0",           
                sascmd="C:/Program Files/SASHome/x86/SASFoundation/9.3/sas.exe") 


# Export data and save as a table file.
write.table(tb1,"tb1.csv",sep=",")

# To A Tab Delimited Text File
write.table(tb1, "tb1.txt", sep="\t") 

# To SAS
library(foreign)
write.foreign(tb1, "tb2.txt", "tb1.sas", package="SAS")

# tb2.txt is the data file
# tb1.sas file is the file you could run in SAS to import this data file in SAS. 
# you can also save ".RData" file and Plot/Graph files with 
# save.image(), save.plot(),
#----------------------------- end import/export -------------------------------

#----------------------------- Functions in R -------------------------------
# R is a high level 4th generation language, i.e. it is a non-procedural, functional language. As such it can be very powerful (if you know how to use it; the next step would be the Matrix...).

#functions in R are used like this:  function_call(arguments)

#many functions are generic to deal with different objects differently, e.g. summary or plot
# just for show, later we will discuss this at length
x <- 1:10                			
y <- 2*x+9+rnorm(10)
linmod <- lm(y~x)

# call the standardized R help file for lm (linear model)
?lm                           
# function_call: summary; argument: data frame df11 output: 
# table for factors and five point summaries for metric variables
summary(df11)    

# function_call: summary; argument: 
# object of type "lm"; same function call, but very different output, a linear model fit
summary(linmod)            		
?summary         #help file for summary

#It is also possible to pass fuctions as arguments to functions, therefore making thing very flexible

x <- rnorm(500)           			#generating 500 realisations of a standard normal distributed random variable 
y <- factor(sample(c(1,2),500,replace=TRUE)) 	#generate 500 labels
spineplot(y~x)           			#make a spineplot, note we have a number of unequal breakpints 
fivenum(x)               			#five point summary or x
spineplot(y~x,breaks=fivenum(x))	 	#we now use the five point summary function as an argument for breakpoints	


# --- Writing functions
#a simple function for calculating mean and variance 

add <- function(a,b) 
  { result = a+b
  return(result)}

v2 <- add(1,2) 
#----------------------------- End functions in R -------------------------------

#----------------------------- Matrix algebra in R ----------------------------
# R can do matrix algebra (nearly as good as matlab)

A <- mat1
B <- t(mat2)+5

A+B             #element wise addition 
A-B             #element wise subtraction 

3*A             #scalar multiplication

A*B             #element wise multiplication, Hadamar product
A%*%t(B)           #matrix multiplication (inner product) 
A%o%B           #outer product, AB'

outer(A,B,FUN="*") #the same; you can use other functions for FUN e.g. sum or "-" 

crossprod(A,B)  #A'B 
crossprod(A)    #A'A

det(crossprod(A))          #determinant of A'A

t(A)            #transpose
diag(5)         #5x5 identity matrix; generic for scalar 
diag(vec1)      #length(vec1) x length(vec1) diagonal matrix with vec1 as diagonal elements; generic for vectors
diag(A)         #returns diagonal elements of A; generic for matrices
sum(diag(A))    #trace of A 

# Hilbert transform
hilbert <- function(n) {i <- 1:n; 1 / outer(i - 1, i, "+") }
(h8 <- hilbert(8))

(sh8 <- solve(h8))

# round and obtain matrix multiply results with 5 decimal places. 
round(sh8%*%h8,digits=5)

solve(crossprod(A))        #Inverse of A if A is square matrix

# Load MASS library
library(MASS)

ginv(crossprod(A))         #Moore-Penrose inverse from package MASS
z <- eigen(crossprod(A))   #eigen structure of matrix
z$val
z$vec

z <- svd(A)     #singular value decomposition of A
z$d             #vector of singular values of A
z$u             #matrix with columns containing left singular vectors of A
z$v             #matrix with columns containing right singular vectors of A

(R <- chol(crossprod(A)))    #Choleski-factorization of A, R'R=A )
z <- qr(A)      #QR decomposition of A
z$qr            #upper triangle that contains decomposition and lower triangle with info on Q decomposition
z$rank          #rank of A
z$qraux         #additional information on Q
z$pivot         #information on pivoting strategy 
# get help about qr
?qr

#convenience functions
rowMeans(A)  #vector of means per row
colMeans(A)  #vector of means per column
rowSums(A)   #vector of sums per row 
colSums(A)   #vector of sums per row

#More specialized functions (BLAS, Lapack, sparse matrices etc.)  in 
#library(Matrix)               

#--------------------- end matrix functions --------------------------------------
#--------------------- GUI in R ---------------------------------------------
# please intall R package
# 'fgui', 

# load library
library(tcltk)
library(fgui)



# simple add function.
add <- function(a,b) 
{ result = a+b
  return(result)}

# obtain results from this function through code
v1 <- add(2,2)
v2 <- add(1,2)
v3 <- add(3,5)

# making GUI for this ADD function
add <- gui(add, argOption=list(CallOption=c("TRUE","FALSE")))


# -------------------------
# some example of GUI by using other two package.
# please intall R package
# 'gWidgetstcltk', 'gWidgets'

# load library
library(gWidgetstcltk)
library(gWidgets)

# Some sample data to test against
x1 <- rnorm(100)
x2 <- runif(100)

# Widget labels
labelX <- "Variable name for data: "
labelY <- "Distribution to compare to: "
labelAlternative <- "One or two sided test?: "
labelP <- "The p-value is: "

# Choices for comboboxes
choicesAlternative <- eval(formals(ks.test)$alternative)
distributions <- c(
  normal = pnorm, 
  exponential = pexp,
  F = pf,
  "log-normal" = plnorm,
  "Student's t" = pt,
  uniform = punif)

# create GUI.
createKsTestGwidgets <- function()
{
  library(gWidgetstcltk)
  options(guiToolkit = "tcltk")
  win <- gwindow("KS Test, gWidgets edition", visible = FALSE)
  
  frmX <- gframe("x", container = win)
  lblX <- glabel(labelX, container = frmX)
  txtX <- gedit(container = frmX)
  
  frmY <- gframe("y", container = win)
  lblY <- glabel(labelY, container = frmY)
  cmbY <- gcombobox(names(distributions), container = frmY)
  
  frmAlternative <- gframe("alternative", container = win)
  lblAlternative <- glabel(labelAlternative, container = frmAlternative)
  cmbAlternative <- gcombobox(choicesAlternative, container = frmAlternative)
  
  btnCalc <- gbutton("Calculate", container = win,
                     handler = function(h, ...)
                     {
                       x <- get(svalue(txtX), mode = "numeric")
                       y <- distributions[[svalue(cmbY)]]
                       alternative <- svalue(cmbAlternative)
                       ans <- ks.test(x, y, alternative = alternative)
                       svalue(txtP) <- format(ans$p.value, digits = 3)
                     }
  )
  frmResults <- gframe("results", container = win)
  lblP <- glabel(labelP, container = frmResults)
  txtP <- gedit(container = frmResults)
  visible(win) <- TRUE
  
}

# launch the function.
createKsTestGwidgets()

#----------------------------- end GUI in R -------------------------------
#---------------------- Fun parts of R ------------------------------------------

#--------------------- Web scraping -------------------------------------- 
## more infos, tutorials, etc.: http://www.omegahat.org/

## Two basic strategies:
## 1) Through API (typically more user-friendly)
## 2) DIY (more data preparation needed)

## Twitter - Getting tweets into R
# install.packages "twitteR"  

# load package
require("twitteR")


# define twitter user
tuser <- getUser("BarackObama")

public.tweets <- publicTimeline()                             ## last 20 public tweets (http://twitter.com/public_timeline)
public.tweets

ob.timeline <- userTimeline("BarackObama", n = 20)            ## last 20 tweets by AW
ob.timeline
## ---------------------------- end webscraping ---------------------------

## ---------------------------- draw a map in R ----------------------------
## install "maps" package and "mapdata" first.

# load package
library(maps)
library(mapdata)

# draw a world map
if (require(maps) && require(mapdata)) {
  par(mar = c(0, 0, 0, 0))
  map("world")
  bd = par()$usr  #get the border
  x = seq(bd[1], bd[2], length = 100)
  y = seq(bd[3], bd[4], length = 100)
  z = sqrt(outer((x - mean(x))^2, (y - mean(y))^2, "+"))
  image(x, y, z, col = cm.colors(25), axes = F, xlab = "", ylab = "")
  map("world", add = TRUE) # add the map to the color image
}

# draw a United States map.
if (require(maps) && require(mapdata)) {
  par(mar = c(0, 0, 0, 0))
  map("state")
  bd = par()$usr  #get the border
  x = seq(bd[1], bd[2], length = 100)
  y = seq(bd[3], bd[4], length = 100)
  z = sqrt(outer((x - mean(x))^2, (y - mean(y))^2, "+"))
  image(x, y, z, col = cm.colors(25), axes = F, xlab = "", ylab = "")
  map("state", add = TRUE) # add the map to the color image
}

## ----------------The fun Package: Use R for Fun! -------------------
## install.packages('fun')

## load "fun" package
library(fun)

if (.Platform$OS.type == "windows") x11() else x11(type = "Xlib")

# Hanoi tower demonstration.
tower_of_hanoi(3)   #  3 layers

## Game 1: Japanese Gomoku game.
gomoku()

## Test Alzheimer's disease by finding out the different character in a character rectangle
x = alzheimer_test()

## Mine Sweeper
if (.Platform$OS.type == "windows") x11() else x11(type = "Xlib")
mine_sweeper()

## ----------------------------- USE R to play some music  -----------------------
#####################################

# Most shining national wind on R
# install library "seewave".

# load library
library(seewave)

fs = 44100; # sample rate 
dt = 1/fs;

T16 = 0.125;

t16=seq(0,T16,by=dt);

k = length(t16);

t4 = seq(0,4*T16,length.out=4*k); 
t8 = seq(0,2*T16,length.out=2*k);

i = length(t4); 
j = length(t8);

# Modification functions 
mod4=(t4^4)*exp(-30*(t4^0.5)); 
mod4=mod4*(1/max(mod4)); 
mod8=(t8^4)*exp(-50*(t8^0.5)); 
mod8=mod8*(1/max(mod8)); 
mod16=(t16^4)*exp(-90*(t16^0.5)); 
mod16=mod16*(1/max(mod16));

f0 = 2*146.8; # reference frequency

ScaleTable = c(2/3,3/4,5/6,15/16,1,9/8,5/4,4/3,3/2,5/3,9/5,15/8,2,9/4,5/2,8/3,3,10/3,15/4,4,1/2,9/16,5/8);

# 1/4 notes 
do0f = mod4*cos(2*pi*ScaleTable[21]*f0*t4); 
re0f = mod4*cos(2*pi*ScaleTable[22]*f0*t4); 
mi0f = mod4*cos(2*pi*ScaleTable[23]*f0*t4);

fa0f = mod4*cos(2*pi*ScaleTable[1]*f0*t4); 
so0f = mod4*cos(2*pi*ScaleTable[2]*f0*t4); 
la0f = mod4*cos(2*pi*ScaleTable[3]*f0*t4); 
ti0f = mod4*cos(2*pi*ScaleTable[4]*f0*t4); 
do1f = mod4*cos(2*pi*ScaleTable[5]*f0*t4); 
re1f = mod4*cos(2*pi*ScaleTable[6]*f0*t4); 
mi1f = mod4*cos(2*pi*ScaleTable[7]*f0*t4); 
fa1f = mod4*cos(2*pi*ScaleTable[8]*f0*t4); 
so1f = mod4*cos(2*pi*ScaleTable[9]*f0*t4); 
la1f = mod4*cos(2*pi*ScaleTable[10]*f0*t4); 
tb1f = mod4*cos(2*pi*ScaleTable[11]*f0*t4); 
ti1f = mod4*cos(2*pi*ScaleTable[12]*f0*t4); 
do2f = mod4*cos(2*pi*ScaleTable[13]*f0*t4); 
re2f = mod4*cos(2*pi*ScaleTable[14]*f0*t4); 
mi2f = mod4*cos(2*pi*ScaleTable[15]*f0*t4); 
fa2f = mod4*cos(2*pi*ScaleTable[16]*f0*t4); 
so2f = mod4*cos(2*pi*ScaleTable[17]*f0*t4); 
la2f = mod4*cos(2*pi*ScaleTable[18]*f0*t4); 
ti2f = mod4*cos(2*pi*ScaleTable[19]*f0*t4); 
do3f = mod4*cos(2*pi*ScaleTable[20]*f0*t4); 
blkf = c (rep(0,i));

# 1/8 notes 
do0e = mod8*cos(2*pi*ScaleTable[21]*f0*t8); 
re0e = mod8*cos(2*pi*ScaleTable[22]*f0*t8); 
mi0e = mod8*cos(2*pi*ScaleTable[23]*f0*t8);

fa0e = mod8*cos(2*pi*ScaleTable[1]*f0*t8); 
so0e = mod8*cos(2*pi*ScaleTable[2]*f0*t8); 
la0e = mod8*cos(2*pi*ScaleTable[3]*f0*t8); 
ti0e = mod8*cos(2*pi*ScaleTable[4]*f0*t8); 
do1e = mod8*cos(2*pi*ScaleTable[5]*f0*t8); 
re1e = mod8*cos(2*pi*ScaleTable[6]*f0*t8); 
mi1e = mod8*cos(2*pi*ScaleTable[7]*f0*t8); 
fa1e = mod8*cos(2*pi*ScaleTable[8]*f0*t8); 
so1e = mod8*cos(2*pi*ScaleTable[9]*f0*t8); 
la1e = mod8*cos(2*pi*ScaleTable[10]*f0*t8); 
tb1e = mod8*cos(2*pi*ScaleTable[11]*f0*t8); 
ti1e = mod8*cos(2*pi*ScaleTable[12]*f0*t8); 
do2e = mod8*cos(2*pi*ScaleTable[13]*f0*t8); 
re2e = mod8*cos(2*pi*ScaleTable[14]*f0*t8); 
mi2e = mod8*cos(2*pi*ScaleTable[15]*f0*t8); 
fa2e = mod8*cos(2*pi*ScaleTable[16]*f0*t8); 
so2e = mod8*cos(2*pi*ScaleTable[17]*f0*t8); 
la2e = mod8*cos(2*pi*ScaleTable[18]*f0*t8); 
ti2e = mod8*cos(2*pi*ScaleTable[19]*f0*t8); 
do3e = mod8*cos(2*pi*ScaleTable[20]*f0*t8); 
blke = c(rep(0,j));

# 1/16 notes 
do0s = mod16*cos(2*pi*ScaleTable[21]*f0*t16); 
re0s = mod16*cos(2*pi*ScaleTable[22]*f0*t16); 
mi0s = mod16*cos(2*pi*ScaleTable[23]*f0*t16);

fa0s = mod16*cos(2*pi*ScaleTable[1]*f0*t16); 
so0s = mod16*cos(2*pi*ScaleTable[2]*f0*t16); 
la0s = mod16*cos(2*pi*ScaleTable[3]*f0*t16); 
ti0s = mod16*cos(2*pi*ScaleTable[4]*f0*t16); 
do1s = mod16*cos(2*pi*ScaleTable[5]*f0*t16); 
re1s = mod16*cos(2*pi*ScaleTable[6]*f0*t16); 
mi1s = mod16*cos(2*pi*ScaleTable[7]*f0*t16); 
fa1s = mod16*cos(2*pi*ScaleTable[8]*f0*t16); 
so1s = mod16*cos(2*pi*ScaleTable[9]*f0*t16); 
la1s = mod16*cos(2*pi*ScaleTable[10]*f0*t16); 
tb1s = mod16*cos(2*pi*ScaleTable[11]*f0*t16); 
ti1s = mod16*cos(2*pi*ScaleTable[12]*f0*t16); 
do2s = mod16*cos(2*pi*ScaleTable[13]*f0*t16); 
re2s = mod16*cos(2*pi*ScaleTable[14]*f0*t16); 
mi2s = mod16*cos(2*pi*ScaleTable[15]*f0*t16); 
fa2s = mod16*cos(2*pi*ScaleTable[16]*f0*t16); 
so2s = mod16*cos(2*pi*ScaleTable[17]*f0*t16); 
la2s = mod16*cos(2*pi*ScaleTable[18]*f0*t16); 
ti2s = mod16*cos(2*pi*ScaleTable[19]*f0*t16); 
do3s = mod16*cos(2*pi*ScaleTable[20]*f0*t16); 
blks = c(rep(0,k));


# Melody by Schau_mal 
part0=c(mi1f,la0e,la0e,do1f,mi1f,
        re1e,re1s,mi1s,re1e,do1e,re1e,do1e,la0f,
        mi1f,la0e,la0e,do1f,mi1f,
        so1e,re1s,mi1s,re1e,do1e,re1e,do1e,ti0e,so0e,
        mi1f,la0e,la0e,do1f,mi1f,
        re1e,re1s,mi1s,re1e,do1e,re1e,do1e,la0e,so0e,
        mi1f,la0e,la0e,do1f,mi1f,
        so1e,mi1e,blkf,blkf,blkf);

part1=c(la0f,la0e,so0e,la0f,la0e,do1e,
        do1f,re1e,do1e,la0f,la0f,
        do1f,do1e,so0e,do1e,re1e,mi1e,so1e,
        so1e,mi1e,re1f,mi1f,mi1f,
        la1e,la1e,la1e,so1e,mi1e,mi1f,do1e,
        la0e,la0e,la0e,mi1e,re1s,mi1s,re1e,re1f,
        mi1e,mi1e,so1e,mi1e,re1e,mi1e,re1e,do1e,
        la0f,so0f,la0f,la0f);

part2=c(mi1e,mi1e,so1e,mi1e,mi1e,so1e,so1e,la1e,
        do2e,la1e,so1f,la1s,do2s,la1e,la1f,
        la0f,la0e,so0e,la0f,do1f,
        re1e,mi1s,re1s,do1e,re1e,mi1f,mi1f,
        la0e,la1e,la1e,so1e,re1e,mi1s,re1s,do1e,re1e,
        mi1f,mi1f,blke,blke,blkf,
        do1e,la0e,la0e,do1e,re1f,so0e,so0e,
        mi1e,so1e,mi1e,re1e,do1f,do1f,
        la0e,do1e,re1e,mi1e,re1e,do1e,so0e,mi0e,
        la0f,la0f,blke,blke,blkf);

part3=c(la0f,la0e,so0e,la0f,do1f,
        re1e,mi1s,re1s,do1e,re1e,mi1f,mi1f,
        la0e,la1e,la1e,so1e,re1e,mi1s,re1s,do1e,re1e,
        mi1f,mi1f,blke,blke,blkf,
        do1e,la0e,la0e,do1e,re1f,so0e,so0e,
        mi1e,so1e,mi1e,re1e,do1f,do1e,do1e,
        la0e,do1e,re1e,mi1e,so1e,mi1e,mi1e,so1e,
        la1f,la1f,la1f,la1f);

part4=c(la1e,la1s,la1s,la1e,la1e,la1e,la1s,so1s,mi1e,re1e,
        re1e,re1s,re1s,mi1e,mi1s,so1s,mi1e,mi1s,re1s,do1e,do1s,la0s,
        la0f,la0e,so0e,la0f,la0e,do1e,
        re1e,mi1s,re1s,do1e,re1e,mi1f,mi1f,
        la1e,so1e,mi1e,re1e,so1e,mi1e,re1e,do1e,
        do1f,do1f,la0s,do1s,re1s,mi1s,re1s,do1s,la0s,do1s);

part5=c(do2e,do2s,do2s,la1e,la1s,la1s,so1e,so1s,so1s,mi1e,mi1s,mi1s,
        re1e,mi1s,re1s,do1e,la0s,so0s,la0s,so0s,do1s,re1s,mi1s,so1s,la1s,re2s,
        do2f,do2f,blks,blks,blks,blks,do1e,re1e,
        mi1f,mi1f,mi1f,so1e,mi1e,
        la1f,la1f,la1e,do1e,so1e,mi1e,
        re1f,re1e,re1s,re1s,re1e,re1e,do1e,re1e,
        mi1f,mi1e,mi1s,mi1s,mi1e,re1s,do1s,ti0e,do1s,re1s,
        mi1f,mi1f,mi1f,so1e,mi1e,
        do2f,la1f,la1f,la1e,do1e,
        re1f,so1f,so1f,la1f,
        ti1f,ti1f,ti1f,ti1f);

part6=c(blkf,blkf,mi1e,so1e,mi1e,so1e,
        mi1f,la0e,la0s,la0s,do1f,la0e,mi1s,la0s,
        do1e,do1s,do1s,re1e,do1s,re1s,mi1f,mi1f,
        mi1f,la0e,la0s,la0s,so1f,re1e,re1s,re1s,
        mi1f,mi1f,mi1s,re1s,do1s,la0s,mi0s,re0s,mi0s,so0s,
        do1f,la0e,la0s,la0s,re1f,so0e,so0s,so0s,
        mi0f,so0e,so0s,so0s,do1f,do1f,
        la0f,do1e,do1s,la0s,mi1e,mi1s,mi1s,re1e,re1s,mi1s);

# Combination, v1 is complete version, v2 is simple version. 
v1=c(part0,part1,part1,part2,part3,part4,part0,part1,part1,part2,part3,part5,part3,part6,part3);
v2=c(part0,part1,part1,part2,part3,part5,part3,part6,part3);

# Let's rock ^_^ 
s = v1; 
s = s/max(s);

listen(s,f=fs)
#------------------------End of Playing Music ---------------------

