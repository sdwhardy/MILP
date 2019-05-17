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

####################################
function lof_dists(ocean,base)
    km=1
    lnth_1deg_LT=111*km
    location[i].coord.x=lof_deg2lgth(location[i].gps.lng-base.lng,lof_lg1deg(location[i].gps.lat+base.lng)*km)
        end
function lof_rotateOcn(ocean)
    theta=atan((ocean.pccs[length(ocean.pccs)].coord.x-ocean.cnces[length(ocean.cnces)].coord.x)/(ocean.cnces[length(ocean.cnces)].coord.y-ocean.pccs[length(ocean.pccs)].coord.y))
    os=0.0
    os=lof_rot(ocean.cnces,theta,os)
    os=lof_rot(ocean.pccs,theta,os)
    return os
end
#loops through to apply appropriate rotations for coords
function lof_rot(location,theta,os)
    for i=1:length(location)
        xy=lof_rotation(location[i].coord.x,location[i].coord.y,theta)
        location[i].coord.x=xy[1]
        location[i].coord.y=xy[2]
        if location[i].coord.x<os
            os=location[i].coord.x
        end
    end
    return os
end
function lof_rotation(x,y,theta)
    co_od=[x y]
    rotated=co_od*[cos(theta) -1*sin(theta);sin(theta) cos(theta)]
    return rotated
end
#translates the entire region
function lof_slideOcn(ocean,os)
    for cnce in ocean.cnces
        cnce.coord.x=lof_shift(cnce.coord.x,abs(os))
        cnce.id="2"*string(cnce.num)
    end
    for pcc in ocean.pccs
        pcc.coord.x=lof_shift(pcc.coord.x,abs(os))
        pcc.id="1"*string(pcc.num)
    end
end
