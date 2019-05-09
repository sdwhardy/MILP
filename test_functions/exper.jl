#=using Gadfly
using DataFrames
df = DataFrame(X = 1:4, Y = [4, 5, 6, 7],C=["A","B","C","D"])
df2 = DataFrame(X = 2:5, Y = [4, 5, 6, 7],C=["E","F","G","H"])

Gadfly.with_theme(:dark) do
    p = layer(df, x=:X, y=:Y, label=:C, Geom.point, Geom.label,Theme(point_label_color="white"))
    p2= layer(df2, x=:X, y=:Y, label=:C, Geom.point, Geom.label,Theme(point_label_color="white"))
    plot(p,p2)
end=#
#Gadfly.pop_theme()
#print(iris)

y = rand(10)
plot(y,annotations=(3,y[3],text("this is #3",:left)),leg=false)
#=annotate!([(5,y[5],text("this is #5",16,:red,:center)),
          (10,y[10],text("this is #10",:right,20,"courier"))])
scatter!(range(2, stop=8, length=6),rand(6),marker=(50,0.2,:orange),
         series_annotations=["series","annotations","map","to","series",
                             text("data",:green)])=#


function cmb(n,k)
 return (factorial(n)/(factorial(k)*factorial(n-k)))
end
function tc(n)
    tot=0.0
    for k=2:n
        tot=tot+cmb(n,k)*k
    end
    return tot
end
println(tc(10))
