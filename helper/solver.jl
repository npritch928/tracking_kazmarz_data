using LinearAlgebra

function gower_Richtarik!(x,A,b,S,Binv)
    StA = BLAS.gemm('T','N',1.,S,A)
    BinvAtS = BLAS.gemm('N','T',1.,Binv,StA)
    Stb = BLAS.gemv('T',1.,S,b)
    #Overwrite S'b with the residual since S'b is no longer needed
    BLAS.gemv!('N',1.,StA,x,-1.,Stb)
    #Create matrix to be used in the pseudoinverse
    cent = BLAS.gemm('N','N',1.,StA,BinvAtS)
    Pseud = pinv(cent)
    Modi = BLAS.gemm('N','N',1.,BinvAtS,Pseud)
    #Update x
    BLAS.gemv!('N',-1.,Modi,Stb,1.,x)
    #Return the residual
    return Stb
end

# Defines the function to calculate the moving average
function move_av(Values::Vector, i::Int64, width::Int64)
    if i > width
        # Average only over necessary width
        av = sum(Values[(i - width + 1): i]) / width
    else
        # Moving average is normal average when below specify width
        av = sum(Values[1:i]) / i
    end
    return av
end
