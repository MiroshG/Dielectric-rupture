using Plots, Random, Statistics, LaTeXStrings, GLM, DataFrames


"################################# NEIGBOURS FUNCTION #################################"

function Dielectric_neighbours(zero_array, potential_matrix) #zero_array contains the elements with zero value in the lattice [i,j] (first row i, second row j)

    neighbours_array=zeros(3,0)  # array with the position of the neighbours (first two columns) and value of the potential (third column)
    
    for k in 1:size(zero_array,2) # do the loop 4 times for each 0 value of the matrix

        i=Int(zero_array[1,k])          # row position of the element with zero potential           zero_array=|i_1, i_2, i_3, ...|
        j=Int(zero_array[2,k])          # column position of the element with zero potential                   |j_1, j_2, j_3, ...|

        # take only the neighbour elements with potential non zero or one 
        if potential_matrix[i-1, j]!= 0.0 && potential_matrix[i-1, j]!= 1.0

            neighbours_array=hcat(neighbours_array, [i-1, j, potential_matrix[i-1, j]])
        end

        if potential_matrix[i+1, j]!= 0.0 && potential_matrix[i+1, j]!= 1.0

            neighbours_array=hcat(neighbours_array, [i+1, j, potential_matrix[i+1, j]])
        end

        if potential_matrix[i, j-1]!= 0.0 && potential_matrix[i, j-1]!= 1.0

            neighbours_array=hcat(neighbours_array, [i, j-1, potential_matrix[i, j-1]])
        end

        if potential_matrix[i, j+1]!= 0.0 && potential_matrix[i, j+1]!= 1.0

            neighbours_array=hcat(neighbours_array, [i, j+1, potential_matrix[i, j+1]])
        end
    end

    return neighbours_array

end

"################################# ZERO POINTS FUNCTION #################################"

function zero_points(potential_matrix)

    zero_array=zeros(2,0)
    

    for i in 1:size(potential_matrix, 1)
        for j in 1:size(potential_matrix, 2)
            if potential_matrix[i,j]==0
                zero_array=hcat(zero_array, [i,j])
            end
            
        end
    end
    return zero_array      # zero_array=|i_1, i_2, i_3, ...|
end                        #            |j_1, j_2, j_3, ...|



"################################# PLOTING BLACK AND WHITE FUNCTION #################################"

function Plotting_black_and_white(potential_matrix)

    plotting_matrix=potential_matrix
    
    for i in 1:size(plotting_matrix,1)
        for j in 1:size(plotting_matrix,2)

            if plotting_matrix[i,j]==0.0 #|| plotting_matrix[i,j]==1.0
                plotting_matrix[i,j]=0.0

            elseif  plotting_matrix[i,j]==1.0
                plotting_matrix[i,j]=0.5

            else 
                plotting_matrix[i,j]=1.0
            end
        end
    end

    return(plotting_matrix)
end



"################################# BOUNDARY 1: POTENTIAL FUNCTION 1 RADIOUS #################################"

function Dielectric_potential_1(L, potential_matrix) # this function calculates the laplace equation for the whole lattice for the 1st boundary condition
    
    # Potential in the center
    center=(L-1)/2+1
    

    "Boundary condition 1 ############################"
    # Potential in the boundary equal to 1
    for i in 1:L
        for j in 1:L
            if (i-center)^2+(j-center)^2>=((L-1)/2)^2
                potential_matrix[i,j]=1
            end
        end
    end

    center=Int(center)
    potential_matrix[center,center]=0

    "#################################################"

    # potential in the rest of the lattice
    previous_potential=similar(potential_matrix)    # different matrix than the potential matrix but with the same size

    while previous_potential != potential_matrix    # loop until the potential reamins constant 

        previous_potential .=potential_matrix       # update the previous potential with the new one

        """for i in 2:L-1
            for j in 2:L-1
                if potential_matrix[i,j]==0         # points with potential equal to zero must remain unchanged
                    potential_matrix[i,j]=0
                elseif potential_matrix[i,j]==1     # points with potential equal to one must remain unchanged
                    potential_matrix[i,j]=1
                else
                    # changes in the potential 
                    potential_matrix[i,j]=1/4*(previous_potential[i+1,j]+previous_potential[i-1,j]+previous_potential[i,j+1]+previous_potential[i,j-1]) 
                end
            end
        end"""
        for i in 2:L-1
            for j in 2:L-1
                if potential_matrix[i,j]==0         # points with potential equal to zero must remain unchanged
                    potential_matrix[i,j]=0
                elseif potential_matrix[i,j]==1     # points with potential equal to one must remain unchanged
                    potential_matrix[i,j]=1
                else
                    # changes in the potential 
                    potential_matrix[i,j]=1/4*(previous_potential[i+1,j]+previous_potential[i-1,j]+previous_potential[i,j+1]+previous_potential[i,j-1]) 
                end
            end
        end
    end

    return(potential_matrix)
end




function Pattern_growth(L,t, eta)  # function to simulate the growth of the dielectric breakdown

    potential_matrix=rand(0.1:1,L,L)  # First lattice to start


    # LOOP FOR THE FIRST BOUNDARY CONDITIONS #############  
    for i in 1:t-1
        
        Dielectric_potential_1(L, potential_matrix) # Equilibrated potential

        # positions of the elements with zero potential (first row are the row values i and second row are the column values j)
        zero_array=zero_points(potential_matrix)  

        # positions of the neighbours elements and their value                                   | i_1, i_2, i_3, ...|
        neighbours_array=Dielectric_neighbours(zero_array, potential_matrix)  # neighbours_array=| j_1, j_2, j_3, ...|       number of columns= 
        # select a random neighbour element                                                      | Φ_1, Φ_2, Φ_3, ...|       number of neighbours
        completed=false

        while (!completed)

            # random neighbour selected (we select a random column)
            selected=Int(rand(1:size(neighbours_array,2)))  
            

            # potential of the selected element divided by the sum of all posible elements
            prob_acceptance=(neighbours_array[3,selected]^eta/sum(neighbours_array[3,:].^eta))   

            if rand()< prob_acceptance
                i=Int(neighbours_array[1, selected])   
                j=Int(neighbours_array[2, selected])
                potential_matrix[i,j]=0
                completed=true
            end

            #for i in 1:size(neighbours_array,2)


        end
        
    end

    return(potential_matrix)
    
end

function Fractal_dimensionality_1( # function to compute the fractal dimensionality. It's defined as the relation between the total number of points "N(r)"
    L,                             # inside a radius "r" and this radius where "D" is de fractal dimensionality as:   
    t, eta)                        #                                  N(r)~r^D

    potential_matrix=Pattern_growth(L,t, eta)

    radius= [x for x in 5:1:exp(3)]   # array with different values of the radius to obtain d with a linear regression # exp(2.5)
    total_length_of_branches=[]        # array to keep the total number of points

    center=(L-1)/2+1             # element of the center of the lattice
    central_r=ones(L,L)

    for R in radius              # loop over the different radius

        for i in 1:L
            for j in 1:L
                if (i-center)^2+(j-center)^2<=R^2
                    central_r[i,j]=potential_matrix[i,j]            #matrix with the central elements of the potential with a radius R
                end
            end
        end

        total_length=size(zero_points(central_r),2) # the number of columns of the zero_points function is the number of elements with zero potential

        push!(total_length_of_branches, total_length)
    end

    total_length_of_branches=convert(Vector{Float64}, total_length_of_branches)
    
    model = lm(@formula(y ~ x), DataFrame(x=log.(radius) , y=log.(total_length_of_branches)))

    fit_params = coef(model)
    slope, intercept = fit_params[2], fit_params[1]

    # Obtain the errors 
    errors = stderror(model)
    slope_error, intercept_error = errors[2], errors[1]

    r_squared=r2(model)

    println("D($eta)= ", slope, " ± ", slope_error )
    println(" ")
    println("r2_D= " , r_squared)
    println(" ")

    return(potential_matrix,total_length_of_branches,radius)

end
################################################################################################################
# POTENTIAL MATRIX AT THE BEGGINING ############################################################################
################################################################################################################

L=101


"""potential_matrix=rand(0.1:1,L,L)

potential_matrix_time_zero= Dielectric_potential_1(L, potential_matrix) 

heatmap(potential_matrix_time_zero, c=:thermal, aspect_ratio=:equal, colorbar=true, xlabel="", ylabel="")
savefig("Color_map_time_zero.png")"""

#################################################################################################################
# PATERN GROWTH #################################################################################################
#################################################################################################################

"""potential_matrix_pattern=Pattern_growth(L,1000, 1.)

heatmap(potential_matrix_pattern, c=:thermal, aspect_ratio=:equal, colorbar=true, xlabel="", ylabel="")

plotting_matrix=Plotting_black_and_white(potential_matrix_pattern)
savefig("Color_map_time_2000.png")

heatmap(plotting_matrix, color=:grays, axis=false, colorbar=false, aspect_ratio=:equal) 

savefig("Dielectric_breakdown_at_t=100.png")"""
    
#################################################################################################################
# FRACTAL DIMENSIONALITY #########################################################################################
#################################################################################################################

eta_array=[0., 0.5, 2., 1.]
time_array=[100, 200, 500, 1000, 2000]
L=101

"""for time in time_array
    eta=0.0

    println("time: ",time)
    potential_matrix_pattern,total_length_of_branches,radius=Fractal_dimensionality_1(L, time, eta)

    plot!(log.(radius), log.(total_length_of_branches), xlabel="ln(r)", ylabel="ln(N)", grid=false, label="time")
    #ylims!(3, 3.5)
    savefig("Evolution_r_N-base_eta.png")

end"""

for eta in eta_array
    if eta!= 2.0
        potential_matrix_pattern,total_length_of_branches,radius=Fractal_dimensionality_1(L, 2000, eta)
    else
        potential_matrix_pattern,total_length_of_branches,radius=Fractal_dimensionality_1(L, 1000, eta)
    end

    plotting_matrix=Plotting_black_and_white(potential_matrix_pattern)

    heatmap(plotting_matrix, color=:grays, axis=false, colorbar=false, aspect_ratio=:equal) 

    savefig("Dielectric_breakdown_$eta.png")

    heatmap(potential_matrix_pattern, c=:thermal, aspect_ratio=:equal, colorbar=true, xlabel="", ylabel="")
    savefig("Color_map_thermal_$eta.png")
end
