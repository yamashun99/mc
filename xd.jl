module Xd
    function _f(x)#被積分関数
        a = 0
        for i in x
         a+=i^2
        end
        return a
    end

    _all_perm(x, n) = Iterators.product([x for i = 1:n]...)

    function Flatt(N,d)#台形公式
        F=0
        for x in _all_perm([i/N for i=0:N-1], d)
            for delta in _all_perm([i/N for i=0:1], d)
                F+=_f(collect(x)+collect(delta))
            end
        end
        F/=N^d*2^d
        return F
    end

    function Fmc(N,d)#モンテカルロ
        F=0
        for ix=1:N^d
            x=[rand() for i=1:d]
            F+=_f(x)
        end
        F/=N^d
        return F
    end

    function sigma(N,d,M)#モンテカルロの分散
        F2=0
        F1=0
        for iM=1:M
            F=Fmc(N,d)
            F2+=F^2
            F1+=F
        end
        F2/=M
        F1/=M
        return sqrt(F2-F1^2)
    end

    function precision(Ns,d,M)
        deltalatt = [(abs(Flatt(N,d)-d/3)/(d/3)) for N in Ns]#台形公式の精度
        deltamc=[(Xd.Fmc(N,d)-d/3)/(d/3) for N in Ns]#モンテカルロの精度
        sigmas = [Xd.sigma(i,d,M) for i in Ns]#モンテカルロの分散
        return deltalatt, deltamc, sigmas
    end
end