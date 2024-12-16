dsem.tt <- list()
#########################################################
dsem.tt[[4]] <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1+l7*y4t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
eta1t2 =~ ly1*y1t2 + ly2*y2t2 + ly3*y3t2+l7*y4t2
eta2t2 =~ ly4*y4t2 + ly5*y5t2 + ly6*y6t2
#
eta1t3 =~ ly1*y1t3 + ly2*y2t3 + ly3*y3t3+l7*y4t3
eta2t3 =~ ly4*y4t3 + ly5*y5t3 + ly6*y6t3
#
eta1t4 =~ ly1*y1t4 + ly2*y2t4 + ly3*y3t4+l7*y4t4
eta2t4 =~ ly4*y4t4 + ly5*y5t4 + ly6*y6t4
#
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
y1t3 ~~ td1*y1t3
y1t4 ~~ td1*y1t4
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
y2t3 ~~ td2*y2t3
y2t4 ~~ td2*y2t4
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
y3t3 ~~ td3*y3t3
y3t4 ~~ td3*y3t4
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
y4t3 ~~ td4*y4t3
y4t4 ~~ td4*y4t4
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
y5t3 ~~ td5*y5t3
y5t4 ~~ td5*y5t4
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
y6t3 ~~ td6*y6t3
y6t4 ~~ td6*y6t4
# 
eta1t2 ~ beta1*eta1t1
eta1t3 ~ beta1*eta1t2
eta1t4 ~ beta1*eta1t3
# 
eta2t2 ~ beta2*eta2t1
eta2t3 ~ beta2*eta2t2
eta2t4 ~ beta2*eta2t3
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
eta1t3 ~~ phi21*eta2t3
eta1t4 ~~ phi21*eta2t4
#
eta1t1 ~~ phi11*eta1t1
eta1t2 ~~ phi11*eta1t2
eta1t3 ~~ phi11*eta1t3
eta1t4 ~~ phi11*eta1t4
#
eta2t1 ~~ phi22*eta2t1
eta2t2 ~~ phi22*eta2t2
eta2t3 ~~ phi22*eta2t3
eta2t4 ~~ phi22*eta2t4
#
u01 =~ lu01*eta1t1+lu01*eta1t2+lu01*eta1t3+lu01*eta1t4
u02 =~ lu02*eta2t1+lu02*eta2t2+lu02*eta2t3+lu02*eta2t4
#
u01 ~~ u02
u01 ~~ u01
u02 ~~ u02
'

#########################################################
dsem.tt[[3]] <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1+l7*y4t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
eta1t2 =~ ly1*y1t2 + ly2*y2t2 + ly3*y3t2+l7*y4t2
eta2t2 =~ ly4*y4t2 + ly5*y5t2 + ly6*y6t2
#
eta1t3 =~ ly1*y1t3 + ly2*y2t3 + ly3*y3t3+l7*y4t3
eta2t3 =~ ly4*y4t3 + ly5*y5t3 + ly6*y6t3
#
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
y1t3 ~~ td1*y1t3
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
y2t3 ~~ td2*y2t3
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
y3t3 ~~ td3*y3t3
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
y4t3 ~~ td4*y4t3
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
y5t3 ~~ td5*y5t3
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
y6t3 ~~ td6*y6t3
# 
eta1t2 ~ beta1*eta1t1
eta1t3 ~ beta1*eta1t2
# 
eta2t2 ~ beta2*eta2t1
eta2t3 ~ beta2*eta2t2
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
eta1t3 ~~ phi21*eta2t3
#
eta1t1 ~~ phi11*eta1t1
eta1t2 ~~ phi11*eta1t2
eta1t3 ~~ phi11*eta1t3
#
eta2t1 ~~ phi22*eta2t1
eta2t2 ~~ phi22*eta2t2
eta2t3 ~~ phi22*eta2t3
#
u01 =~ lu01*eta1t1+lu01*eta1t2+lu01*eta1t3
u02 =~ lu02*eta2t1+lu02*eta2t2+lu02*eta2t3
#
u01 ~~ u02
u01 ~~ u01
u02 ~~ u02
'

#########################################################
dsem.tt[[2]] <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1+l7*y4t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
eta1t2 =~ ly1*y1t2 + ly2*y2t2 + ly3*y3t2+l7*y4t2
eta2t2 =~ ly4*y4t2 + ly5*y5t2 + ly6*y6t2
#
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
# 
eta1t2 ~ beta1*eta1t1
# 
eta2t2 ~ beta2*eta2t1
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
#
eta1t1 ~~ phi11*eta1t1
eta1t2 ~~ phi11*eta1t2
#
eta2t1 ~~ phi22*eta2t1
eta2t2 ~~ phi22*eta2t2
#
u01 =~ lu01*eta1t1+lu01*eta1t2
u02 =~ lu02*eta2t1+lu02*eta2t2
#
u01 ~~ u02
u01 ~~ u01
u02 ~~ u02
'


#########################################################
dsem.tt[[1]] <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1+l7*y4t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
y1t1 ~~ td1*y1t1
#
y2t1 ~~ td2*y2t1
#
y3t1 ~~ td3*y3t1
#
y4t1 ~~ td4*y4t1
#
y5t1 ~~ td5*y5t1
#
y6t1 ~~ td6*y6t1
# 
eta1t1 ~~ phi21*eta2t1
'
#########################################################
dsem.tt[[5]] <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1+l7*y4t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
eta1t2 =~ ly1*y1t2 + ly2*y2t2 + ly3*y3t2+l7*y4t2
eta2t2 =~ ly4*y4t2 + ly5*y5t2 + ly6*y6t2
#
eta1t3 =~ ly1*y1t3 + ly2*y2t3 + ly3*y3t3+l7*y4t3
eta2t3 =~ ly4*y4t3 + ly5*y5t3 + ly6*y6t3
#
eta1t4 =~ ly1*y1t4 + ly2*y2t4 + ly3*y3t4+l7*y4t4
eta2t4 =~ ly4*y4t4 + ly5*y5t4 + ly6*y6t4
#
eta1t5 =~ ly1*y1t5 + ly2*y2t5 + ly3*y3t5+l7*y4t5
eta2t5 =~ ly4*y4t5 + ly5*y5t5 + ly6*y6t5
#
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
y1t3 ~~ td1*y1t3
y1t4 ~~ td1*y1t4
y1t5 ~~ td1*y1t5
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
y2t3 ~~ td2*y2t3
y2t4 ~~ td2*y2t4
y2t5 ~~ td2*y2t5
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
y3t3 ~~ td3*y3t3
y3t4 ~~ td3*y3t4
y3t5 ~~ td3*y3t5
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
y4t3 ~~ td4*y4t3
y4t4 ~~ td4*y4t4
y4t5 ~~ td4*y4t5
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
y5t3 ~~ td5*y5t3
y5t4 ~~ td5*y5t4
y5t5 ~~ td5*y5t5
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
y6t3 ~~ td6*y6t3
y6t4 ~~ td6*y6t4
y6t5 ~~ td6*y6t5
# 
eta1t2 ~ beta1*eta1t1
eta1t3 ~ beta1*eta1t2
eta1t4 ~ beta1*eta1t3
eta1t5 ~ beta1*eta1t4
# 
eta2t2 ~ beta2*eta2t1
eta2t3 ~ beta2*eta2t2
eta2t4 ~ beta2*eta2t3
eta2t5 ~ beta2*eta2t4
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
eta1t3 ~~ phi21*eta2t3
eta1t4 ~~ phi21*eta2t4
eta1t5 ~~ phi21*eta2t5
#
eta1t1 ~~ phi11*eta1t1
eta1t2 ~~ phi11*eta1t2
eta1t3 ~~ phi11*eta1t3
eta1t4 ~~ phi11*eta1t4
eta1t5 ~~ phi11*eta1t5
#
eta2t1 ~~ phi22*eta2t1
eta2t2 ~~ phi22*eta2t2
eta2t3 ~~ phi22*eta2t3
eta2t4 ~~ phi22*eta2t4
eta2t5 ~~ phi22*eta2t5
#
u01 =~ lu01*eta1t1+lu01*eta1t2+lu01*eta1t3+lu01*eta1t4+lu01*eta1t5
u02 =~ lu02*eta2t1+lu02*eta2t2+lu02*eta2t3+lu02*eta2t4+lu02*eta2t5
#
u01 ~~ u02
u01 ~~ u01
u02 ~~ u02
'

#########################################################
dsem.tt[[10]] <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1+l7*y4t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
eta1t2 =~ ly1*y1t2 + ly2*y2t2 + ly3*y3t2+l7*y4t2
eta2t2 =~ ly4*y4t2 + ly5*y5t2 + ly6*y6t2
#
eta1t3 =~ ly1*y1t3 + ly2*y2t3 + ly3*y3t3+l7*y4t3
eta2t3 =~ ly4*y4t3 + ly5*y5t3 + ly6*y6t3
#
eta1t4 =~ ly1*y1t4 + ly2*y2t4 + ly3*y3t4+l7*y4t4
eta2t4 =~ ly4*y4t4 + ly5*y5t4 + ly6*y6t4
#
eta1t5 =~ ly1*y1t5 + ly2*y2t5 + ly3*y3t5+l7*y4t5
eta2t5 =~ ly4*y4t5 + ly5*y5t5 + ly6*y6t5
#
eta1t6 =~ ly1*y1t6 + ly2*y2t6 + ly3*y3t6+l7*y4t6
eta2t6 =~ ly4*y4t6 + ly5*y5t6 + ly6*y6t6
#
eta1t7 =~ ly1*y1t7 + ly2*y2t7 + ly3*y3t7+l7*y4t7
eta2t7 =~ ly4*y4t7 + ly5*y5t7 + ly6*y6t7
#
eta1t8 =~ ly1*y1t8 + ly2*y2t8 + ly3*y3t8+l7*y4t8
eta2t8 =~ ly4*y4t8 + ly5*y5t8 + ly6*y6t8
#
eta1t9 =~ ly1*y1t9 + ly2*y2t9 + ly3*y3t9+l7*y4t9
eta2t9 =~ ly4*y4t9 + ly5*y5t9 + ly6*y6t9
#
eta1t10 =~ ly1*y1t10 + ly2*y2t10 + ly3*y3t10+l7*y4t10
eta2t10 =~ ly4*y4t10 + ly5*y5t10 + ly6*y6t10
#
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
y1t3 ~~ td1*y1t3
y1t4 ~~ td1*y1t4
y1t5 ~~ td1*y1t5
y1t6 ~~ td1*y1t6
y1t7 ~~ td1*y1t7
y1t8 ~~ td1*y1t8
y1t9 ~~ td1*y1t9
y1t10 ~~ td1*y1t10
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
y2t3 ~~ td2*y2t3
y2t4 ~~ td2*y2t4
y2t5 ~~ td2*y2t5
y2t6 ~~ td2*y2t6
y2t7 ~~ td2*y2t7
y2t8 ~~ td2*y2t8
y2t9 ~~ td2*y2t9
y2t10 ~~ td2*y2t10
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
y3t3 ~~ td3*y3t3
y3t4 ~~ td3*y3t4
y3t5 ~~ td3*y3t5
y3t6 ~~ td3*y3t6
y3t7 ~~ td3*y3t7
y3t8 ~~ td3*y3t8
y3t9 ~~ td3*y3t9
y3t10 ~~ td3*y3t10
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
y4t3 ~~ td4*y4t3
y4t4 ~~ td4*y4t4
y4t5 ~~ td4*y4t5
y4t6 ~~ td4*y4t6
y4t7 ~~ td4*y4t7
y4t8 ~~ td4*y4t8
y4t9 ~~ td4*y4t9
y4t10 ~~ td4*y4t10
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
y5t3 ~~ td5*y5t3
y5t4 ~~ td5*y5t4
y5t5 ~~ td5*y5t5
y5t6 ~~ td5*y5t6
y5t7 ~~ td5*y5t7
y5t8 ~~ td5*y5t8
y5t9 ~~ td5*y5t9
y5t10 ~~ td5*y5t10
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
y6t3 ~~ td6*y6t3
y6t4 ~~ td6*y6t4
y6t5 ~~ td6*y6t5
y6t6 ~~ td6*y6t6
y6t7 ~~ td6*y6t7
y6t8 ~~ td6*y6t8
y6t9 ~~ td6*y6t9
y6t10 ~~ td6*y6t10
# 
eta1t2 ~ beta1*eta1t1
eta1t3 ~ beta1*eta1t2
eta1t4 ~ beta1*eta1t3
eta1t5 ~ beta1*eta1t4
eta1t6 ~ beta1*eta1t5
eta1t7 ~ beta1*eta1t6
eta1t8 ~ beta1*eta1t7
eta1t9 ~ beta1*eta1t8
eta1t10 ~ beta1*eta1t9
# 
eta2t2 ~ beta2*eta2t1
eta2t3 ~ beta2*eta2t2
eta2t4 ~ beta2*eta2t3
eta2t5 ~ beta2*eta2t4
eta2t6 ~ beta2*eta2t5
eta2t7 ~ beta2*eta2t6
eta2t8 ~ beta2*eta2t7
eta2t9 ~ beta2*eta2t8
eta2t10 ~ beta2*eta2t9
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
eta1t3 ~~ phi21*eta2t3
eta1t4 ~~ phi21*eta2t4
eta1t5 ~~ phi21*eta2t5
eta1t6 ~~ phi21*eta2t6
eta1t7 ~~ phi21*eta2t7
eta1t8 ~~ phi21*eta2t8
eta1t9 ~~ phi21*eta2t9
eta1t10 ~~ phi21*eta2t10
#
eta1t1 ~~ phi11*eta1t1
eta1t2 ~~ phi11*eta1t2
eta1t3 ~~ phi11*eta1t3
eta1t4 ~~ phi11*eta1t4
eta1t5 ~~ phi11*eta1t5
eta1t6 ~~ phi11*eta1t6
eta1t7 ~~ phi11*eta1t7
eta1t8 ~~ phi11*eta1t8
eta1t9 ~~ phi11*eta1t9
eta1t10 ~~ phi11*eta1t10
#
eta2t1 ~~ phi22*eta2t1
eta2t2 ~~ phi22*eta2t2
eta2t3 ~~ phi22*eta2t3
eta2t4 ~~ phi22*eta2t4
eta2t5 ~~ phi22*eta2t5
eta2t6 ~~ phi22*eta2t6
eta2t7 ~~ phi22*eta2t7
eta2t8 ~~ phi22*eta2t8
eta2t9 ~~ phi22*eta2t9
eta2t10 ~~ phi22*eta2t10
#
u01 =~ lu01*eta1t1+lu01*eta1t2+lu01*eta1t3+lu01*eta1t4+lu01*eta1t5+lu01*eta1t6+lu01*eta1t7+lu01*eta1t8+lu01*eta1t9+lu01*eta1t10
u02 =~ lu02*eta2t1+lu02*eta2t2+lu02*eta2t3+lu02*eta2t4+lu02*eta2t5+lu02*eta2t6+lu02*eta2t7+lu02*eta2t8+lu02*eta2t9+lu02*eta2t10
#
u01 ~~ u02
u01 ~~ u01
u02 ~~ u02
'

#########################################################
dsem.tt[[15]] <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1+l7*y4t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
eta1t2 =~ ly1*y1t2 + ly2*y2t2 + ly3*y3t2+l7*y4t2
eta2t2 =~ ly4*y4t2 + ly5*y5t2 + ly6*y6t2
#
eta1t3 =~ ly1*y1t3 + ly2*y2t3 + ly3*y3t3+l7*y4t3
eta2t3 =~ ly4*y4t3 + ly5*y5t3 + ly6*y6t3
#
eta1t4 =~ ly1*y1t4 + ly2*y2t4 + ly3*y3t4+l7*y4t4
eta2t4 =~ ly4*y4t4 + ly5*y5t4 + ly6*y6t4
#
eta1t5 =~ ly1*y1t5 + ly2*y2t5 + ly3*y3t5+l7*y4t5
eta2t5 =~ ly4*y4t5 + ly5*y5t5 + ly6*y6t5
#
eta1t6 =~ ly1*y1t6 + ly2*y2t6 + ly3*y3t6+l7*y4t6
eta2t6 =~ ly4*y4t6 + ly5*y5t6 + ly6*y6t6
#
eta1t7 =~ ly1*y1t7 + ly2*y2t7 + ly3*y3t7+l7*y4t7
eta2t7 =~ ly4*y4t7 + ly5*y5t7 + ly6*y6t7
#
eta1t8 =~ ly1*y1t8 + ly2*y2t8 + ly3*y3t8+l7*y4t8
eta2t8 =~ ly4*y4t8 + ly5*y5t8 + ly6*y6t8
#
eta1t9 =~ ly1*y1t9 + ly2*y2t9 + ly3*y3t9+l7*y4t9
eta2t9 =~ ly4*y4t9 + ly5*y5t9 + ly6*y6t9
#
eta1t10 =~ ly1*y1t10 + ly2*y2t10 + ly3*y3t10+l7*y4t10
eta2t10 =~ ly4*y4t10 + ly5*y5t10 + ly6*y6t10
#
eta1t11 =~ ly1*y1t11 + ly2*y2t11 + ly3*y3t11+l7*y4t11
eta2t11 =~ ly4*y4t11 + ly5*y5t11 + ly6*y6t11
#
eta1t12 =~ ly1*y1t12 + ly2*y2t12 + ly3*y3t12+l7*y4t12
eta2t12 =~ ly4*y4t12 + ly5*y5t12 + ly6*y6t12
#
eta1t13 =~ ly1*y1t13 + ly2*y2t13 + ly3*y3t13+l7*y4t13
eta2t13 =~ ly4*y4t13 + ly5*y5t13 + ly6*y6t13
#
eta1t14 =~ ly1*y1t14 + ly2*y2t14 + ly3*y3t14+l7*y4t14
eta2t14 =~ ly4*y4t14 + ly5*y5t14 + ly6*y6t14
#
eta1t15 =~ ly1*y1t15 + ly2*y2t15 + ly3*y3t15+l7*y4t15
eta2t15 =~ ly4*y4t15 + ly5*y5t15 + ly6*y6t15
#
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
y1t3 ~~ td1*y1t3
y1t4 ~~ td1*y1t4
y1t5 ~~ td1*y1t5
y1t6 ~~ td1*y1t6
y1t7 ~~ td1*y1t7
y1t8 ~~ td1*y1t8
y1t9 ~~ td1*y1t9
y1t10 ~~ td1*y1t10
y1t11 ~~ td1*y1t11
y1t12 ~~ td1*y1t12
y1t13 ~~ td1*y1t13
y1t14 ~~ td1*y1t14
y1t15 ~~ td1*y1t15
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
y2t3 ~~ td2*y2t3
y2t4 ~~ td2*y2t4
y2t5 ~~ td2*y2t5
y2t6 ~~ td2*y2t6
y2t7 ~~ td2*y2t7
y2t8 ~~ td2*y2t8
y2t9 ~~ td2*y2t9
y2t10 ~~ td2*y2t10
y2t11 ~~ td2*y2t11
y2t12 ~~ td2*y2t12
y2t13 ~~ td2*y2t13
y2t14 ~~ td2*y2t14
y2t15 ~~ td2*y2t15
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
y3t3 ~~ td3*y3t3
y3t4 ~~ td3*y3t4
y3t5 ~~ td3*y3t5
y3t6 ~~ td3*y3t6
y3t7 ~~ td3*y3t7
y3t8 ~~ td3*y3t8
y3t9 ~~ td3*y3t9
y3t10 ~~ td3*y3t10
y3t11 ~~ td3*y3t11
y3t12 ~~ td3*y3t12
y3t13 ~~ td3*y3t13
y3t14 ~~ td3*y3t14
y3t15 ~~ td3*y3t15
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
y4t3 ~~ td4*y4t3
y4t4 ~~ td4*y4t4
y4t5 ~~ td4*y4t5
y4t6 ~~ td4*y4t6
y4t7 ~~ td4*y4t7
y4t8 ~~ td4*y4t8
y4t9 ~~ td4*y4t9
y4t10 ~~ td4*y4t10
y4t11 ~~ td4*y4t11
y4t12 ~~ td4*y4t12
y4t13 ~~ td4*y4t13
y4t14 ~~ td4*y4t14
y4t15 ~~ td4*y4t15
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
y5t3 ~~ td5*y5t3
y5t4 ~~ td5*y5t4
y5t5 ~~ td5*y5t5
y5t6 ~~ td5*y5t6
y5t7 ~~ td5*y5t7
y5t8 ~~ td5*y5t8
y5t9 ~~ td5*y5t9
y5t10 ~~ td5*y5t10
y5t11 ~~ td5*y5t11
y5t12 ~~ td5*y5t12
y5t13 ~~ td5*y5t13
y5t14 ~~ td5*y5t14
y5t15 ~~ td5*y5t15
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
y6t3 ~~ td6*y6t3
y6t4 ~~ td6*y6t4
y6t5 ~~ td6*y6t5
y6t6 ~~ td6*y6t6
y6t7 ~~ td6*y6t7
y6t8 ~~ td6*y6t8
y6t9 ~~ td6*y6t9
y6t10 ~~ td6*y6t10
y6t11 ~~ td6*y6t11
y6t12 ~~ td6*y6t12
y6t13 ~~ td6*y6t13
y6t14 ~~ td6*y6t14
y6t15 ~~ td6*y6t15
# 
eta1t2 ~ beta1*eta1t1
eta1t3 ~ beta1*eta1t2
eta1t4 ~ beta1*eta1t3
eta1t5 ~ beta1*eta1t4
eta1t6 ~ beta1*eta1t5
eta1t7 ~ beta1*eta1t6
eta1t8 ~ beta1*eta1t7
eta1t9 ~ beta1*eta1t8
eta1t10 ~ beta1*eta1t9
eta1t11 ~ beta1*eta1t10
eta1t12 ~ beta1*eta1t11
eta1t13 ~ beta1*eta1t12
eta1t14 ~ beta1*eta1t13
eta1t15 ~ beta1*eta1t14
# 
eta2t2 ~ beta2*eta2t1
eta2t3 ~ beta2*eta2t2
eta2t4 ~ beta2*eta2t3
eta2t5 ~ beta2*eta2t4
eta2t6 ~ beta2*eta2t5
eta2t7 ~ beta2*eta2t6
eta2t8 ~ beta2*eta2t7
eta2t9 ~ beta2*eta2t8
eta2t10 ~ beta2*eta2t9
eta2t11 ~ beta2*eta2t10
eta2t12 ~ beta2*eta2t11
eta2t13 ~ beta2*eta2t12
eta2t14 ~ beta2*eta2t13
eta2t15 ~ beta2*eta2t14
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
eta1t3 ~~ phi21*eta2t3
eta1t4 ~~ phi21*eta2t4
eta1t5 ~~ phi21*eta2t5
eta1t6 ~~ phi21*eta2t6
eta1t7 ~~ phi21*eta2t7
eta1t8 ~~ phi21*eta2t8
eta1t9 ~~ phi21*eta2t9
eta1t10 ~~ phi21*eta2t10
eta1t11 ~~ phi21*eta2t11
eta1t12 ~~ phi21*eta2t12
eta1t13 ~~ phi21*eta2t13
eta1t14 ~~ phi21*eta2t14
eta1t15 ~~ phi21*eta2t15
#
eta1t1 ~~ phi11*eta1t1
eta1t2 ~~ phi11*eta1t2
eta1t3 ~~ phi11*eta1t3
eta1t4 ~~ phi11*eta1t4
eta1t5 ~~ phi11*eta1t5
eta1t6 ~~ phi11*eta1t6
eta1t7 ~~ phi11*eta1t7
eta1t8 ~~ phi11*eta1t8
eta1t9 ~~ phi11*eta1t9
eta1t10 ~~ phi11*eta1t10
eta1t11 ~~ phi11*eta1t11
eta1t12 ~~ phi11*eta1t12
eta1t13 ~~ phi11*eta1t13
eta1t14 ~~ phi11*eta1t14
eta1t15 ~~ phi11*eta1t15
#
eta2t1 ~~ phi22*eta2t1
eta2t2 ~~ phi22*eta2t2
eta2t3 ~~ phi22*eta2t3
eta2t4 ~~ phi22*eta2t4
eta2t5 ~~ phi22*eta2t5
eta2t6 ~~ phi22*eta2t6
eta2t7 ~~ phi22*eta2t7
eta2t8 ~~ phi22*eta2t8
eta2t9 ~~ phi22*eta2t9
eta2t10 ~~ phi22*eta2t10
eta2t11 ~~ phi22*eta2t11
eta2t12 ~~ phi22*eta2t12
eta2t13 ~~ phi22*eta2t13
eta2t14 ~~ phi22*eta2t14
eta2t15 ~~ phi22*eta2t15
#
u01 =~ lu01*eta1t1+lu01*eta1t2+lu01*eta1t3+lu01*eta1t4+lu01*eta1t5+lu01*eta1t6+lu01*eta1t7+lu01*eta1t8+lu01*eta1t9+lu01*eta1t10+
       lu01*eta1t11+lu01*eta1t12+lu01*eta1t13+lu01*eta1t14+lu01*eta1t15
u02 =~ lu02*eta2t1+lu02*eta2t2+lu02*eta2t3+lu02*eta2t4+lu02*eta2t5+lu02*eta2t6+lu02*eta2t7+lu02*eta2t8+lu02*eta2t9+lu02*eta2t10+
       lu02*eta2t11+lu02*eta2t12+lu02*eta2t13+lu02*eta2t14+lu02*eta2t15
        
#
u01 ~~ u02
u01 ~~ u01
u02 ~~ u02
'


#########################################################
dsem.tt[[30]] <- '
eta1t1 =~ ly1*y1t1 + ly2*y2t1 + ly3*y3t1+l7*y4t1
eta2t1 =~ ly4*y4t1 + ly5*y5t1 + ly6*y6t1
#
eta1t2 =~ ly1*y1t2 + ly2*y2t2 + ly3*y3t2+l7*y4t2
eta2t2 =~ ly4*y4t2 + ly5*y5t2 + ly6*y6t2
#
eta1t3 =~ ly1*y1t3 + ly2*y2t3 + ly3*y3t3+l7*y4t3
eta2t3 =~ ly4*y4t3 + ly5*y5t3 + ly6*y6t3
#
eta1t4 =~ ly1*y1t4 + ly2*y2t4 + ly3*y3t4+l7*y4t4
eta2t4 =~ ly4*y4t4 + ly5*y5t4 + ly6*y6t4
#
eta1t5 =~ ly1*y1t5 + ly2*y2t5 + ly3*y3t5+l7*y4t5
eta2t5 =~ ly4*y4t5 + ly5*y5t5 + ly6*y6t5
#
eta1t6 =~ ly1*y1t6 + ly2*y2t6 + ly3*y3t6+l7*y4t6
eta2t6 =~ ly4*y4t6 + ly5*y5t6 + ly6*y6t6
#
eta1t7 =~ ly1*y1t7 + ly2*y2t7 + ly3*y3t7+l7*y4t7
eta2t7 =~ ly4*y4t7 + ly5*y5t7 + ly6*y6t7
#
eta1t8 =~ ly1*y1t8 + ly2*y2t8 + ly3*y3t8+l7*y4t8
eta2t8 =~ ly4*y4t8 + ly5*y5t8 + ly6*y6t8
#
eta1t9 =~ ly1*y1t9 + ly2*y2t9 + ly3*y3t9+l7*y4t9
eta2t9 =~ ly4*y4t9 + ly5*y5t9 + ly6*y6t9
#
eta1t10 =~ ly1*y1t10 + ly2*y2t10 + ly3*y3t10+l7*y4t10
eta2t10 =~ ly4*y4t10 + ly5*y5t10 + ly6*y6t10
#
eta1t11 =~ ly1*y1t11 + ly2*y2t11 + ly3*y3t11+l7*y4t11
eta2t11 =~ ly4*y4t11 + ly5*y5t11 + ly6*y6t11
#
eta1t12 =~ ly1*y1t12 + ly2*y2t12 + ly3*y3t12+l7*y4t12
eta2t12 =~ ly4*y4t12 + ly5*y5t12 + ly6*y6t12
#
eta1t13 =~ ly1*y1t13 + ly2*y2t13 + ly3*y3t13+l7*y4t13
eta2t13 =~ ly4*y4t13 + ly5*y5t13 + ly6*y6t13
#
eta1t14 =~ ly1*y1t14 + ly2*y2t14 + ly3*y3t14+l7*y4t14
eta2t14 =~ ly4*y4t14 + ly5*y5t14 + ly6*y6t14
#
eta1t15 =~ ly1*y1t15 + ly2*y2t15 + ly3*y3t15+l7*y4t15
eta2t15 =~ ly4*y4t15 + ly5*y5t15 + ly6*y6t15
#
eta1t16 =~ ly1*y1t16 + ly2*y2t16 + ly3*y3t16+l7*y4t16
eta2t16 =~ ly4*y4t16 + ly5*y5t16 + ly6*y6t16
#
eta1t17 =~ ly1*y1t17 + ly2*y2t17 + ly3*y3t17+l7*y4t17
eta2t17 =~ ly4*y4t17 + ly5*y5t17 + ly6*y6t17
#
eta1t18 =~ ly1*y1t18 + ly2*y2t18 + ly3*y3t18+l7*y4t18
eta2t18 =~ ly4*y4t18 + ly5*y5t18 + ly6*y6t18
#
eta1t19 =~ ly1*y1t19 + ly2*y2t19 + ly3*y3t19+l7*y4t19
eta2t19 =~ ly4*y4t19 + ly5*y5t19 + ly6*y6t19
#
eta1t20 =~ ly1*y1t20 + ly2*y2t20 + ly3*y3t20+l7*y4t20
eta2t20 =~ ly4*y4t20 + ly5*y5t20 + ly6*y6t20
#
eta1t21 =~ ly1*y1t21 + ly2*y2t21 + ly3*y3t21+l7*y4t21
eta2t21 =~ ly4*y4t21 + ly5*y5t21 + ly6*y6t21
#
eta1t22 =~ ly1*y1t22 + ly2*y2t22 + ly3*y3t22+l7*y4t22
eta2t22 =~ ly4*y4t22 + ly5*y5t22 + ly6*y6t22
#
eta1t23 =~ ly1*y1t23 + ly2*y2t23 + ly3*y3t23+l7*y4t23
eta2t23 =~ ly4*y4t23 + ly5*y5t23 + ly6*y6t23
#
eta1t24 =~ ly1*y1t24 + ly2*y2t24 + ly3*y3t24+l7*y4t24
eta2t24 =~ ly4*y4t24 + ly5*y5t24 + ly6*y6t24
#
eta1t25 =~ ly1*y1t25 + ly2*y2t25 + ly3*y3t25+l7*y4t25
eta2t25 =~ ly4*y4t25 + ly5*y5t25 + ly6*y6t25
#
eta1t26 =~ ly1*y1t26 + ly2*y2t26 + ly3*y3t26+l7*y4t26
eta2t26 =~ ly4*y4t26 + ly5*y5t26 + ly6*y6t26
#
eta1t27 =~ ly1*y1t27 + ly2*y2t27 + ly3*y3t27+l7*y4t27
eta2t27 =~ ly4*y4t27 + ly5*y5t27 + ly6*y6t27
#
eta1t28 =~ ly1*y1t28 + ly2*y2t28 + ly3*y3t28+l7*y4t28
eta2t28 =~ ly4*y4t28 + ly5*y5t28 + ly6*y6t28
#
eta1t29 =~ ly1*y1t29 + ly2*y2t29 + ly3*y3t29+l7*y4t29
eta2t29 =~ ly4*y4t29 + ly5*y5t29 + ly6*y6t29
#
eta1t30 =~ ly1*y1t30 + ly2*y2t30 + ly3*y3t30+l7*y4t30
eta2t30 =~ ly4*y4t30 + ly5*y5t30 + ly6*y6t30
#
y1t1 ~~ td1*y1t1
y1t2 ~~ td1*y1t2
y1t3 ~~ td1*y1t3
y1t4 ~~ td1*y1t4
y1t5 ~~ td1*y1t5
y1t6 ~~ td1*y1t6
y1t7 ~~ td1*y1t7
y1t8 ~~ td1*y1t8
y1t9 ~~ td1*y1t9
y1t10 ~~ td1*y1t10
y1t11 ~~ td1*y1t11
y1t12 ~~ td1*y1t12
y1t13 ~~ td1*y1t13
y1t14 ~~ td1*y1t14
y1t15 ~~ td1*y1t15
y1t16 ~~ td1*y1t16
y1t17 ~~ td1*y1t17
y1t18 ~~ td1*y1t18
y1t19 ~~ td1*y1t19
y1t20 ~~ td1*y1t20
y1t21 ~~ td1*y1t21
y1t22 ~~ td1*y1t22
y1t23 ~~ td1*y1t23
y1t24 ~~ td1*y1t24
y1t25 ~~ td1*y1t25
y1t26 ~~ td1*y1t26
y1t27 ~~ td1*y1t27
y1t28 ~~ td1*y1t28
y1t29 ~~ td1*y1t29
y1t30 ~~ td1*y1t30
#
y2t1 ~~ td2*y2t1
y2t2 ~~ td2*y2t2
y2t3 ~~ td2*y2t3
y2t4 ~~ td2*y2t4
y2t5 ~~ td2*y2t5
y2t6 ~~ td2*y2t6
y2t7 ~~ td2*y2t7
y2t8 ~~ td2*y2t8
y2t9 ~~ td2*y2t9
y2t10 ~~ td2*y2t10
y2t11 ~~ td2*y2t11
y2t12 ~~ td2*y2t12
y2t13 ~~ td2*y2t13
y2t14 ~~ td2*y2t14
y2t15 ~~ td2*y2t15
y2t16 ~~ td2*y2t16
y2t17 ~~ td2*y2t17
y2t18 ~~ td2*y2t18
y2t19 ~~ td2*y2t19
y2t20 ~~ td2*y2t20
y2t21 ~~ td2*y2t21
y2t22 ~~ td2*y2t22
y2t23 ~~ td2*y2t23
y2t24 ~~ td2*y2t24
y2t25 ~~ td2*y2t25
y2t26 ~~ td2*y2t26
y2t27 ~~ td2*y2t27
y2t28 ~~ td2*y2t28
y2t29 ~~ td2*y2t29
y2t30 ~~ td2*y2t30
#
y3t1 ~~ td3*y3t1
y3t2 ~~ td3*y3t2
y3t3 ~~ td3*y3t3
y3t4 ~~ td3*y3t4
y3t5 ~~ td3*y3t5
y3t6 ~~ td3*y3t6
y3t7 ~~ td3*y3t7
y3t8 ~~ td3*y3t8
y3t9 ~~ td3*y3t9
y3t10 ~~ td3*y3t10
y3t11 ~~ td3*y3t11
y3t12 ~~ td3*y3t12
y3t13 ~~ td3*y3t13
y3t14 ~~ td3*y3t14
y3t15 ~~ td3*y3t15
y3t16 ~~ td3*y3t16
y3t17 ~~ td3*y3t17
y3t18 ~~ td3*y3t18
y3t19 ~~ td3*y3t19
y3t20 ~~ td3*y3t20
y3t21 ~~ td3*y3t21
y3t22 ~~ td3*y3t22
y3t23 ~~ td3*y3t23
y3t24 ~~ td3*y3t24
y3t25 ~~ td3*y3t25
y3t26 ~~ td3*y3t26
y3t27 ~~ td3*y3t27
y3t28 ~~ td3*y3t28
y3t29 ~~ td3*y3t29
y3t30 ~~ td3*y3t30
#
y4t1 ~~ td4*y4t1
y4t2 ~~ td4*y4t2
y4t3 ~~ td4*y4t3
y4t4 ~~ td4*y4t4
y4t5 ~~ td4*y4t5
y4t6 ~~ td4*y4t6
y4t7 ~~ td4*y4t7
y4t8 ~~ td4*y4t8
y4t9 ~~ td4*y4t9
y4t10 ~~ td4*y4t10
y4t11 ~~ td4*y4t11
y4t12 ~~ td4*y4t12
y4t13 ~~ td4*y4t13
y4t14 ~~ td4*y4t14
y4t15 ~~ td4*y4t15
y4t16 ~~ td4*y4t16
y4t17 ~~ td4*y4t17
y4t18 ~~ td4*y4t18
y4t19 ~~ td4*y4t19
y4t20 ~~ td4*y4t20
y4t21 ~~ td4*y4t21
y4t22 ~~ td4*y4t22
y4t23 ~~ td4*y4t23
y4t24 ~~ td4*y4t24
y4t25 ~~ td4*y4t25
y4t26 ~~ td4*y4t26
y4t27 ~~ td4*y4t27
y4t28 ~~ td4*y4t28
y4t29 ~~ td4*y4t29
y4t30 ~~ td4*y4t30
#
y5t1 ~~ td5*y5t1
y5t2 ~~ td5*y5t2
y5t3 ~~ td5*y5t3
y5t4 ~~ td5*y5t4
y5t5 ~~ td5*y5t5
y5t6 ~~ td5*y5t6
y5t7 ~~ td5*y5t7
y5t8 ~~ td5*y5t8
y5t9 ~~ td5*y5t9
y5t10 ~~ td5*y5t10
y5t11 ~~ td5*y5t11
y5t12 ~~ td5*y5t12
y5t13 ~~ td5*y5t13
y5t14 ~~ td5*y5t14
y5t15 ~~ td5*y5t15
y5t16 ~~ td5*y5t16
y5t17 ~~ td5*y5t17
y5t18 ~~ td5*y5t18
y5t19 ~~ td5*y5t19
y5t20 ~~ td5*y5t20
y5t21 ~~ td5*y5t21
y5t22 ~~ td5*y5t22
y5t23 ~~ td5*y5t23
y5t24 ~~ td5*y5t24
y5t25 ~~ td5*y5t25
y5t26 ~~ td5*y5t26
y5t27 ~~ td5*y5t27
y5t28 ~~ td5*y5t28
y5t29 ~~ td5*y5t29
y5t30 ~~ td5*y5t30
#
y6t1 ~~ td6*y6t1
y6t2 ~~ td6*y6t2
y6t3 ~~ td6*y6t3
y6t4 ~~ td6*y6t4
y6t5 ~~ td6*y6t5
y6t6 ~~ td6*y6t6
y6t7 ~~ td6*y6t7
y6t8 ~~ td6*y6t8
y6t9 ~~ td6*y6t9
y6t10 ~~ td6*y6t10
y6t11 ~~ td6*y6t11
y6t12 ~~ td6*y6t12
y6t13 ~~ td6*y6t13
y6t14 ~~ td6*y6t14
y6t15 ~~ td6*y6t15
y6t16 ~~ td6*y6t16
y6t17 ~~ td6*y6t17
y6t18 ~~ td6*y6t18
y6t19 ~~ td6*y6t19
y6t20 ~~ td6*y6t20
y6t21 ~~ td6*y6t21
y6t22 ~~ td6*y6t22
y6t23 ~~ td6*y6t23
y6t24 ~~ td6*y6t24
y6t25 ~~ td6*y6t25
y6t26 ~~ td6*y6t26
y6t27 ~~ td6*y6t27
y6t28 ~~ td6*y6t28
y6t29 ~~ td6*y6t29
y6t30 ~~ td6*y6t30
# 
eta1t2 ~ beta1*eta1t1
eta1t3 ~ beta1*eta1t2
eta1t4 ~ beta1*eta1t3
eta1t5 ~ beta1*eta1t4
eta1t6 ~ beta1*eta1t5
eta1t7 ~ beta1*eta1t6
eta1t8 ~ beta1*eta1t7
eta1t9 ~ beta1*eta1t8
eta1t10 ~ beta1*eta1t9
eta1t11 ~ beta1*eta1t10
eta1t12 ~ beta1*eta1t11
eta1t13 ~ beta1*eta1t12
eta1t14 ~ beta1*eta1t13
eta1t15 ~ beta1*eta1t14
eta1t16 ~ beta1*eta1t15
eta1t17 ~ beta1*eta1t16
eta1t18 ~ beta1*eta1t17
eta1t19 ~ beta1*eta1t18
eta1t20 ~ beta1*eta1t19
eta1t21 ~ beta1*eta1t20
eta1t22 ~ beta1*eta1t21
eta1t23 ~ beta1*eta1t22
eta1t24 ~ beta1*eta1t23
eta1t25 ~ beta1*eta1t24
eta1t26 ~ beta1*eta1t25
eta1t27 ~ beta1*eta1t26
eta1t28 ~ beta1*eta1t27
eta1t29 ~ beta1*eta1t28
eta1t30 ~ beta1*eta1t29
# 
eta2t2 ~ beta2*eta2t1
eta2t3 ~ beta2*eta2t2
eta2t4 ~ beta2*eta2t3
eta2t5 ~ beta2*eta2t4
eta2t6 ~ beta2*eta2t5
eta2t7 ~ beta2*eta2t6
eta2t8 ~ beta2*eta2t7
eta2t9 ~ beta2*eta2t8
eta2t10 ~ beta2*eta2t9
eta2t11 ~ beta2*eta2t10
eta2t12 ~ beta2*eta2t11
eta2t13 ~ beta2*eta2t12
eta2t14 ~ beta2*eta2t13
eta2t15 ~ beta2*eta2t14
eta2t16 ~ beta2*eta2t15
eta2t17 ~ beta2*eta2t16
eta2t18 ~ beta2*eta2t17
eta2t19 ~ beta2*eta2t18
eta2t20 ~ beta2*eta2t19
eta2t21 ~ beta2*eta2t20
eta2t22 ~ beta2*eta2t21
eta2t23 ~ beta2*eta2t22
eta2t24 ~ beta2*eta2t23
eta2t25 ~ beta2*eta2t24
eta2t26 ~ beta2*eta2t25
eta2t27 ~ beta2*eta2t26
eta2t28 ~ beta2*eta2t27
eta2t29 ~ beta2*eta2t28
eta2t30 ~ beta2*eta2t29
# 
eta1t1 ~~ phi21*eta2t1
eta1t2 ~~ phi21*eta2t2
eta1t3 ~~ phi21*eta2t3
eta1t4 ~~ phi21*eta2t4
eta1t5 ~~ phi21*eta2t5
eta1t6 ~~ phi21*eta2t6
eta1t7 ~~ phi21*eta2t7
eta1t8 ~~ phi21*eta2t8
eta1t9 ~~ phi21*eta2t9
eta1t10 ~~ phi21*eta2t10
eta1t11 ~~ phi21*eta2t11
eta1t12 ~~ phi21*eta2t12
eta1t13 ~~ phi21*eta2t13
eta1t14 ~~ phi21*eta2t14
eta1t15 ~~ phi21*eta2t15
eta1t16 ~~ phi21*eta2t16
eta1t17 ~~ phi21*eta2t17
eta1t18 ~~ phi21*eta2t18
eta1t19 ~~ phi21*eta2t19
eta1t20 ~~ phi21*eta2t20
eta1t21 ~~ phi21*eta2t21
eta1t22 ~~ phi21*eta2t22
eta1t23 ~~ phi21*eta2t23
eta1t24 ~~ phi21*eta2t24
eta1t25 ~~ phi21*eta2t25
eta1t26 ~~ phi21*eta2t26
eta1t27 ~~ phi21*eta2t27
eta1t28 ~~ phi21*eta2t28
eta1t29 ~~ phi21*eta2t29
eta1t30 ~~ phi21*eta2t30
#
eta1t1 ~~ phi11*eta1t1
eta1t2 ~~ phi11*eta1t2
eta1t3 ~~ phi11*eta1t3
eta1t4 ~~ phi11*eta1t4
eta1t5 ~~ phi11*eta1t5
eta1t6 ~~ phi11*eta1t6
eta1t7 ~~ phi11*eta1t7
eta1t8 ~~ phi11*eta1t8
eta1t9 ~~ phi11*eta1t9
eta1t10 ~~ phi11*eta1t10
eta1t11 ~~ phi11*eta1t11
eta1t12 ~~ phi11*eta1t12
eta1t13 ~~ phi11*eta1t13
eta1t14 ~~ phi11*eta1t14
eta1t15 ~~ phi11*eta1t15
eta1t16 ~~ phi11*eta1t16
eta1t17 ~~ phi11*eta1t17
eta1t18 ~~ phi11*eta1t18
eta1t19 ~~ phi11*eta1t19
eta1t20 ~~ phi11*eta1t20
eta1t21 ~~ phi11*eta1t21
eta1t22 ~~ phi11*eta1t22
eta1t23 ~~ phi11*eta1t23
eta1t24 ~~ phi11*eta1t24
eta1t25 ~~ phi11*eta1t25
eta1t26 ~~ phi11*eta1t26
eta1t27 ~~ phi11*eta1t27
eta1t28 ~~ phi11*eta1t28
eta1t29 ~~ phi11*eta1t29
eta1t30 ~~ phi11*eta1t30
#
eta2t1 ~~ phi22*eta2t1
eta2t2 ~~ phi22*eta2t2
eta2t3 ~~ phi22*eta2t3
eta2t4 ~~ phi22*eta2t4
eta2t5 ~~ phi22*eta2t5
eta2t6 ~~ phi22*eta2t6
eta2t7 ~~ phi22*eta2t7
eta2t8 ~~ phi22*eta2t8
eta2t9 ~~ phi22*eta2t9
eta2t10 ~~ phi22*eta2t10
eta2t11 ~~ phi22*eta2t11
eta2t12 ~~ phi22*eta2t12
eta2t13 ~~ phi22*eta2t13
eta2t14 ~~ phi22*eta2t14
eta2t15 ~~ phi22*eta2t15
eta2t16 ~~ phi22*eta2t16
eta2t17 ~~ phi22*eta2t17
eta2t18 ~~ phi22*eta2t18
eta2t19 ~~ phi22*eta2t19
eta2t20 ~~ phi22*eta2t20
eta2t21 ~~ phi22*eta2t21
eta2t22 ~~ phi22*eta2t22
eta2t23 ~~ phi22*eta2t23
eta2t24 ~~ phi22*eta2t24
eta2t25 ~~ phi22*eta2t25
eta2t26 ~~ phi22*eta2t26
eta2t27 ~~ phi22*eta2t27
eta2t28 ~~ phi22*eta2t28
eta2t29 ~~ phi22*eta2t29
eta2t30 ~~ phi22*eta2t30
#
u01 =~ lu01*eta1t1+lu01*eta1t2+lu01*eta1t3+lu01*eta1t4+lu01*eta1t5+lu01*eta1t6+lu01*eta1t7+lu01*eta1t8+lu01*eta1t9+lu01*eta1t10+
       lu01*eta1t11+lu01*eta1t12+lu01*eta1t13+lu01*eta1t14+lu01*eta1t15+lu01*eta1t16+lu01*eta1t17+lu01*eta1t18+lu01*eta1t19+lu01*eta1t20+
       lu01*eta1t21+lu01*eta1t22+lu01*eta1t23+lu01*eta1t24+lu01*eta1t25+lu01*eta1t26+lu01*eta1t27+lu01*eta1t28+lu01*eta1t29+lu01*eta1t30
u02 =~ lu02*eta2t1+lu02*eta2t2+lu02*eta2t3+lu02*eta2t4+lu02*eta2t5+lu02*eta2t6+lu02*eta2t7+lu02*eta2t8+lu02*eta2t9+lu02*eta2t10+
       lu02*eta2t11+lu02*eta2t12+lu02*eta2t13+lu02*eta2t14+lu02*eta2t15+lu02*eta2t16+lu02*eta2t17+lu02*eta2t18+lu02*eta2t19+lu02*eta2t20+
       lu02*eta2t21+lu02*eta2t22+lu02*eta2t23+lu02*eta2t24+lu02*eta2t25+lu02*eta2t26+lu02*eta2t27+lu02*eta2t28+lu02*eta2t29+lu02*eta2t30
        
#
u01 ~~ u02
u01 ~~ u01
u02 ~~ u02
'



