from scipy.spatial import distance
import math
import random
import copy
import time as tm

def euclidean_dist(point1, point2):
    '''
    Returns the rounded Euclidean distance between two points.

    Parameters
    ----------
    point1 : tuple
        The first point.
    point2 : tuple
        The second point.
    
    Returns
    -------
    float
        The Euclidean distance between the two points.
    '''
    return math.ceil(distance.euclidean(point1,point2))

def eucl(point1,point2):
    '''
    Returns not rounded Euclidean distance between two points.
    '''
    return distance.euclidean(point1,point2)

class Customer:
    '''
    Class representing a customer.

    Attributes
    ----------
    customer_id : int
        The ID of the customer.
    x_coord : int
        The x coordinate of the customer.
    y_coord : int   
        The y coordinate of the customer.
    demand : int
        The demand of the customer.
    ready_time : int
        The ready time of the customer.
    due_time : int
        The due time of the customer.
    service_time : int
        The service time of the customer.
    distance_to_depot : int
        The distance between the customer and the depot.
    
    Methods
    -------
    distance(customer)
        Returns the distance between the customer and another customer.
    __repr__()
        Returns the string representation of the customer.
    '''
    def __init__(self, data):
        self.customer_id = data[0]
        self.x_coord = int(data[1])
        self.y_coord = int(data[2])
        self.demand = int(data[3])
        self.ready_time = int(data[4])
        self.due_time = int(data[5])
        self.service_time = int(data[6])
        self.distance_to_depot = 0

    def distance(self,customer):
        '''
        Returns the distance (rounded) between the customer and another customer.
        '''
        return euclidean_dist((self.x_coord, self.y_coord), (customer.x_coord, customer.y_coord))
    
    def distance_not_rounded(self,customer):
        '''
        Returns the distance between the customer and another customer.
        '''
        return eucl((self.x_coord, self.y_coord), (customer.x_coord, customer.y_coord))
    
    def __repr__(self):
        return "Customer: " + str(self.customer_id)
    
    def __gt__(self, other):
        return self.customer_id > other.customer_id
    
class Vehicle:
    '''
    Class representing a vehicle.

    Attributes
    ----------
    id : int
        The ID of the vehicle.
    customers : list
        The list of customers in the vehicle.
    capacity : int
        The capacity of the vehicle.
    time : int
        The time of the vehicle.
    distance : int
        The distance of the vehicle.
    max_time : int
        The maximum time of the vehicle.
    time_over : int
        The time over of the vehicle.
    
    Methods
    -------
    add(customer, index = -1)
        Adds a customer to the vehicle.
    remove(customer)
        Removes a customer from the vehicle.
    swap_customers(customer1, customer2)
        Swaps two customers in the vehicle.
    recalculate()
        Recalculates the time, distance and time over of the vehicle.
    reconfigure(customer1, customer2)
        Reconfigures the vehicle.
    copy()
        Returns a copy of the vehicle.
    __repr__()
        Returns the string representation of the vehicle.
    '''
    def __init__(self, id, depot):
        self.id = id
        self.customers = [depot,depot]
        self.capacity = 0
        self.time = 0
        self.distance = 0
        self.max_time = depot.due_time
        self.time_over = 0
    
    def add(self, customer, index = -1):
        '''
        Adds a customer to the vehicle.
        
        Parameters
        ----------
        customer : Customer
            The customer to be added.
        index : int
            The index at which the customer should be added.
        
        Returns
        -------
        None
        '''
        self.capacity += customer.demand
        self.customers.insert(index, customer)
        self.recalculate()

    def remove(self, customer):
        '''
        Removes a customer from the vehicle.

        Parameters
        ----------
        customer : Customer
            The customer to be removed.
        
        Returns
        -------
        int
            The index of the removed customer.
        '''
        index = self.customers.index(customer)
        self.customers.remove(customer)
        self.capacity -= customer.demand
        self.recalculate()
        return index
    
    def swap_customers(self, customer1, customer2):
        '''
        Swaps two customers in the vehicle.
        
        Parameters
        ----------
        customer1 : Customer
            The first customer to be swapped.
        customer2 : Customer
            The second customer to be swapped.
        
        Returns
        -------
        None
        '''
        index = self.remove(customer1)
        self.add(customer2, index)
        self.recalculate()

    def recalculate(self):
        '''
        Recalculates the time, distance and time over of the vehicle.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        self.time = 0
        self.distance = 0
        time_over = 0
        for i in range(len(self.customers[1:])):
            curr = self.customers[i]
            previous = self.customers[i-1]
            dist = curr.distance(previous)
            self.distance += dist
            self.time += dist + previous.service_time
            if curr.ready_time > self.time:
                self.time = curr.ready_time
            time_over += max(0, self.time- curr.due_time)
        self.time_over = time_over

    def reconfigure(self,customer1, customer2):
        '''
        Reconfigures the vehicle.

        Parameters
        ----------
        customer1 : Customer
            The first customer to be reconfigured.
        customer2 : Customer
            The second customer to be reconfigured.
        
        Returns
        -------
        None
        '''
        index1 = self.customers.index(customer1)
        index2 = self.customers.index(customer2)
        self.customers[index1] = customer2
        self.customers[index2] = customer1
        self.recalculate()

    def copy(self):
        return copy.deepcopy(self)

    def __repr__(self):
        return "Vehicle: " + str(self.id)
    

class Solution:
    '''
    Class representing a solution.

    Attributes
    ----------
    vehicles : dict
        The dictionary of vehicles in the solution.
    customers : dict
        The dictionary of customers in the solution.
    removed : list
        The list of removed vehicles.
    
    Methods
    -------
    add(vehicle)
        Adds a vehicle to the solution.
    remove(index)
        Removes a vehicle from the solution.
    swap_random_customers()
        Swaps two random customers in the solution.
    move_to_existing_vehicle()
        Moves a customer to an existing vehicle.
    move_to_new_vehicle()
        Moves a customer to a new vehicle.
    get_neighbor()
        Returns a random neighbor of the solution.
    copy()
        Returns a copy of the solution.
    fitness(capacity, num_of_vehicles)
        Returns the fitness of the solution.
    is_feasible(capacity, num_of_vehicles)
        Returns whether the solution is feasible.
    '''
    def __init__(self, data = None) -> None:
        # dictionary vehicle.id:vehicle
        self.vehicles = {}
        # dictionary Customer:vehicle.id
        self.customers = {}
        self.removed = []

    def add(self,vehicle):
        '''
        Adds a vehicle to the solution.
        
        Parameters
        ----------
        vehicle : Vehicle
            The vehicle to be added.
        
        Returns
        -------
        None
        '''
        self.vehicles[vehicle.id] = vehicle
        for customer in vehicle.customers[1:-1]:
            self.customers[customer] = vehicle.id

    def remove(self, index):
        '''
        Removes a vehicle from the solution.

        Parameters
        ----------
        index : int
            The index of the vehicle to be removed.
        
        Returns
        -------
        None
        '''
        self.vehicles.pop(index)

    # options for neighborhood:
    #   - swap 2 customers - LOWER POSSIBILITY
    #       - swap customers between 2 vehicles
    #       - change order of customers in route
    #   - move customer from 1 vehicle to another - HIGHER POSSIBILITY
    #       - remove from 1 route and add to another
    #       - add another vehicle and route
    #       - remove existing vechicle and route
    def swap_random_customers(self):
        '''
        Swaps two random customers in the solution.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        '''
        cust1, cust2 = random.sample(list(self.customers.keys()), 2)
        veh1 = self.customers[cust1]
        veh2 = self.customers[cust2]

        if veh1 != veh2:
            self.vehicles[veh1].swap_customers(cust1, cust2)
            self.vehicles[veh2].swap_customers(cust2, cust1)
            self.customers[cust1] = veh2
            self.customers[cust2] = veh1
        else:
            self.vehicles[veh1].reconfigure(cust1, cust2)
    
    def move_to_existing_vehicle(self):
        '''
        Moves a customer to an existing vehicle.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        '''
        customer = random.sample(list(self.customers.keys()), 1)[0]
        current_vehicle = self.customers[customer]
        new_vehicle = None
        while new_vehicle != current_vehicle:
            new_vehicle = random.sample(list(self.customers.values()),1)[0]
        new_index = random.randrange(1,len(self.vehicles[new_vehicle].customers)-1)
        self.vehicles[new_vehicle].add(customer,new_index)
        self.vehicles[current_vehicle].remove(customer)
        if len(self.vehicles[current_vehicle].customers) == 2:
            self.remove(current_vehicle)
            self.removed.append(current_vehicle.id)
        self.customers[customer] = self.vehicles[new_vehicle].id

    def move_to_new_vehicle(self):
        '''
        Moves a customer to a new vehicle.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        # change how depot is memorised
        customer = random.sample(list(self.customers.keys()), 1)[0]
        current_vehicle = self.customers[customer]
        new_id = 0
        existing_id = list(self.vehicles.keys())[0]
        vehicle = None
        if len(self.removed) == 0:
            new_id = max(list(self.vehicles.keys())) +1
        else:
            new_id = self.removed[0]
            self.removed.remove(new_id)

        vehicle = Vehicle(new_id, self.vehicles[existing_id].customers[0])
        self.vehicles[current_vehicle].customers.remove(customer)
        vehicle.add(customer)
        self.vehicles[new_id] = vehicle
        self.customers[customer] = new_id
        if len(self.vehicles[current_vehicle].customers) == 2:
            del self.vehicles[current_vehicle]
            self.removed.append(current_vehicle)

    def get_neighbor(self):
        '''
        Returns a random neighbor of the solution.

        Parameters
        ----------
        None

        Returns
        -------
        Solution
            The random neighbor of the solution.
        '''
        p = random.random()
        copy = self.copy()
        if p < 0.65:
            if p < 0.1:
                copy.move_to_new_vehicle()
            else:
                copy.move_to_existing_vehicle()
        else:
            copy.swap_random_customers()

        return copy

    def copy(self):
        return copy.deepcopy(self)
    
            
    def fitness(self, capacity, num_of_vehicles):
        '''
        Returns the fitness of the solution.
        
        Parameters
        ----------
        capacity : int
            The capacity of the vehicle.
        num_of_vehicles : int
            The number of vehicles.
        
        Returns
        -------
        list
            The list of parameters of the fitness function.
        '''
        time_over = 0
        capacity_over = 0
        vehicles_over = max(0,len(self.vehicles) - num_of_vehicles)
        sum_of_distances = 0
        for vehicle in self.vehicles.values():
            if vehicle.capacity > capacity:
                capacity_over += vehicle.capacity - capacity
            if vehicle.time_over > 0:
                time_over += vehicle.time_over
            sum_of_distances += vehicle.distance
        return [num_of_vehicles, sum_of_distances, time_over, capacity_over, vehicles_over]
    
    def is_feasible(self, capacity, num_of_vehicles):
        '''
        Returns whether the solution is feasible.
        
        Parameters
        ----------
        capacity : int
            The capacity of the vehicle.
        num_of_vehicles : int
            The number of vehicles.
        
        Returns
        -------
        bool
            Whether the solution is feasible.
        '''
        for vehicle in self.vehicles.values():
            if vehicle.capacity > capacity:
                #print("CAPACITY OVER")
                return False
            elif vehicle.time_over > 0:
                #print("TIME OVER")
                return False
            elif len(self.vehicles) > num_of_vehicles:
                #print("VEHICLES OVER")
                return False
            elif vehicle.time > vehicle.customers[0].due_time:
                #print("DUE TIME OVER")
                return False
            
        return True



class Problem:
    '''
    Class representing a problem.
    
    Attributes
    ----------
    solution_vehicles : int
        The number of vehicles in the solution.
    vehicle_capacity : int
        The capacity of the vehicle.
    customers : list
        The list of customers in the problem.
    importance_factor : int
        The importance factor of the fitness function.
    alpha : int
        The alpha parameter of the fitness function.
    beta : int
        The beta parameter of the fitness function.
    gamma : int
        The gamma parameter of the fitness function.
    penalty_factor : int
        The penalty factor of the fitness function.
    
    Methods
    -------
    fitness(solution)
        Returns the fitness of the solution.
    is_feasible(solution)
        Returns whether the solution is feasible.
    list_representation(solution)
        Returns the list representation of the solution.
    update_params(solution)
        Updates the parameters of the fitness function.
    '''
    def __init__(self, filename):

        self.solution_vehicles = 0
        self.vehicle_capacity = 0
        self.obj_function_calls = 0
        self.importance_factor = 20000
        self.alpha = 20
        self.beta = 20
        self.gamma = 200
        self.penalty_factor = 1.5


        with open(filename) as file:
            rows = file.readlines()
            self.solution_vehicles, self.vehicle_capacity = [int(i) for i in " ".join(rows[2].split()).split(" ")]
            self.customers = []
            for customer in rows[7:]:
                row = " ".join(customer.split()).split(" ")
                if len(row) < 7:
                    break
                self.customers.append(Customer(row))

    def fitness(self, solution):
        '''
        Returns the fitness of the solution.
        
        Parameters
        ----------
        solution : Solution
            The solution.
        
        Returns
        -------
        float
            The fitness of the solution.
        '''
        self.obj_function_calls += 1
        params = solution.fitness(self.vehicle_capacity, self.solution_vehicles)
        return params[0] * self.importance_factor + params[1] + self.alpha * params[2] + self.beta* params[3]  + self.gamma * params[4]
        

    def is_feasible(self,solution):
        '''
        Returns whether the solution is feasible.
        
        Parameters
        ----------
        solution : Solution
            The solution.
        
        Returns
        -------
        bool
            Whether the solution is feasible.
        '''
        return solution.is_feasible(self.vehicle_capacity, self.solution_vehicles)
    
    def list_representation(self,solution):
        '''
        Returns the list representation of the solution.
        
        Parameters
        ----------
        solution : Solution
            The solution.
        
        Returns
        -------
        list
            The list representation of the solution.
        '''
        list_rep = list(solution.vehicles.values())
        list_rep = [x.customers for x in list_rep]
        list_rep.sort(key=lambda x: x[1])
        return list_rep

    def update_params(self,solution):
        '''
        Updates the parameters of the fitness function.
        
        Parameters
        ----------
        solution : Solution
            The solution.
        
        Returns
        -------
        None
        '''
        if self.is_feasible(solution):
            self.alpha = 10
            self.beta = 10
            self.gamma = 200
        else:
            params = solution.fitness(self.vehicle_capacity, self.solution_vehicles)
            if params[2] > 0:
                self.alpha *= self.penalty_factor
            if params[3] > 0:
                self.beta *= self.penalty_factor
            if params[4] > 0:
                self.gamma *= self.penalty_factor
        



if __name__ == '__main__':

    problem = Problem("data/i5.txt")
    # ===============================================
    # GREEDY ALGORITHM
    # ===============================================
    # creating initial solution with greedy approach
    # algorithm creates solution based on distance and due_time
    # 
    solution = Solution()
    depot = problem.customers[0]
    customers = problem.customers[1:]
    for i,customer in enumerate(customers):
        customer.distance_to_depot = customer.distance(depot)
    customers.sort(key=lambda x: x.due_time - x.distance_to_depot)
    id = 0
    while len(customers) != 0:
        vehicle = Vehicle(id, depot)
        time = 0
        for customer in customers.copy():
            time_to_add = vehicle.customers[-2].service_time + vehicle.customers[-2].distance(customer)
            if time + time_to_add < customer.due_time:
                if vehicle.capacity + customer.demand < problem.vehicle_capacity:
                    vehicle.add(customer)
                    customers.remove(customer)
                    if customer.ready_time < time+time_to_add:
                        time += time_to_add
                    else:
                        time = customer.ready_time
                else:
                    break
        solution.add(vehicle)
        id += 1
    # ===============================================
    # TABU SEARCH
    # ===============================================
    time_defined = 300
    exe_time = tm.time() + time_defined
    tabu_tenure = 3
    tabu_iterations = 1000
    tabu_list = []

    best_solution = solution.copy()
    best_fitness = problem.fitness(solution)
    best_num_of_vehicles = len(solution.vehicles)

    neighbourhood_size = 300
    best_neighbours_size = 5
    print(best_fitness)
    counter =0 
    for i in range(tabu_iterations):
        print(i/tabu_iterations)
        if(exe_time < tm.time()) or counter > 50:
            break
        # generating neighbourhood and choosing random from list of best neighbours
        neighbours = []
        flag = False
        starting_best_fitness = best_fitness
        for j in range(neighbourhood_size):
            neighbour = solution.get_neighbor()
            while problem.list_representation(neighbour) in tabu_list:
                neighbour = solution.get_neighbor()
            neighbour_fitness = problem.fitness(neighbour)
            neighbours.append((neighbour, neighbour_fitness))
            if neighbour_fitness < best_fitness and problem.is_feasible(neighbour):
                best_solution = neighbour
                best_fitness = neighbour_fitness
                best_num_of_vehicles = len(neighbour.vehicles)
                best_neighbour = neighbour
                print(best_fitness)
                flag = True
        if flag:
            counter = 0
            solution = neighbour
            if len(tabu_list) == tabu_tenure:
                tabu_list.pop(0)
            tabu_list.append(problem.list_representation(best_neighbour))
            continue

        neighbours.sort(key = lambda x : x[1])


        best_neighbour = neighbours[:best_neighbours_size][random.randrange(0,best_neighbours_size)][0]
        solution = best_neighbour

        if len(tabu_list) == tabu_tenure:
            tabu_list.pop(0)
        tabu_list.append(problem.list_representation(best_neighbour))

        best_neighbour_fitness = problem.fitness(best_neighbour)
        if best_neighbour_fitness < best_fitness and problem.is_feasible(best_neighbour) and len(best_neighbour.vehicles) <= best_num_of_vehicles:
                best_fitness = best_neighbour_fitness
                best_num_of_vehicles = len(best_neighbour.vehicles)
                best_solution = best_neighbour
        if starting_best_fitness == best_fitness:
            counter +=1
        else:
            counter = 0
        print(best_fitness)
    if(time_defined//60 != 1 and time_defined//60 != 5):
        f = open("results/res-"+"un" + "-i5" + ".txt", "w")
    else:
        f = open("results/res-"+str(time_defined//60) + "m" + "-i5" + ".txt", "w")
    id = 0
    f.write(str(len(best_solution.vehicles))+"\n")
    sum_of_distance = 0
    for vehicles in best_solution.vehicles.values():
        time = 0
        customers = vehicles.customers
        f.write(str(id)+": "+str(customers[0].customer_id)+"("+str(time)+")"+"->")
        for i in range(len(customers)-1):
            sum_of_distance += customers[i+1].distance_not_rounded(customers[i])
            time += customers[i+1].distance(customers[i])
            if time < customers[i+1].ready_time:
                time = customers[i+1].ready_time
            f.write(str(customers[i+1].customer_id)+"("+str(round(time, 2))+")")
            if i != len(customers)-2:
                f.write("->")
            time += customers[i+1].service_time
        f.write("\n")
        id += 1
    f.write(str(round(sum_of_distance, 2)))