import pandas as pd
import numpy as np
import copy
import random as rd
from collections import defaultdict
import matplotlib.pyplot as plt
import time


###################################
def check_connection(area1,area2,region,neighbors):
    unvisited = copy.deepcopy(region)
    visited = []
    queue = []
    reach = 0
    
    visited.append(area1)
    unvisited.pop(unvisited.index(area1))
    queue.append(area1)
    
    while queue != []:
#         print("visited: ", visited)
#         print("unvisited: ", unvisited)
#         print("queue: ", queue)
        nei = intersection(neighbors[queue[0]],unvisited)
#         print("node:", queue[0])
#         print("Downstream: ", nei)
        queue.pop(0)
        for aaa in nei:
            visited.append(aaa)
            unvisited.pop(unvisited.index(aaa))
            queue.append(aaa)
            if aaa == area2:
#                 print(area1, "is able to reach", area2, "!!!")
                reach = 1
        if reach == 1:
            break
    return reach



###################################
def min_g_A(N, R_k, d, A):   
    g_min = np.inf  
    for i in N:        
        g = 0       
        for j in R_k:
            g += d[A.index(i),A.index(j)]            
        if g < g_min:
            g_min = g
            min_g_i = i            
#     print(min_g_i)       
    return min_g_i



###################################
def AssignEnclaves(partition,d,A,neighbors):
    parti = copy.deepcopy(partition)
    A_a = []
    for rrr in partition:
        for aaa in rrr:
            A_a.append(aaa)
    e = []
    for aaa in A:
        if aaa not in A_a:
            e.append(aaa)
    e_ori = copy.deepcopy(e)
    
    while e != []:
        A_i = select(e,A_a,neighbors)
        i = A.index(A_i)
        
        ita = region_share_border_area(parti, A_i, neighbors)
      
        R_k = min_g_R(parti, ita, A_i, d, A)
        
        R_k_new = copy.deepcopy(R_k)
        R_k_new.append(A_i)
        
        parti.pop(parti.index(R_k))
        parti.append(R_k_new)
        
        A_a.append(A_i)
                
        e.pop(e.index(A_i))
        
    P_feasible = copy.deepcopy(parti)
    P_feasible = list(np.sort(P_feasible))
    
    return P_feasible, e_ori



###################################
def select(e,A_a,neighbors):
    neib = set()
    for aaa in A_a:
        neib.update(set(neighbors[aaa]))
    to_choose_from = list(neib.intersection(set(e)))    
    A_i = to_choose_from[rd.randrange(0, len(to_choose_from), 1)]
    return A_i



###################################
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3



###################################
def min_g_R(parti, ita, A_i, d,A):    
    g_min = np.inf    
    for rrr in ita:       
        region = parti[rrr]        
        g = 0        
        for j in region:
            if A.index(A_i) <= j:
                g += d[A.index(A_i),A.index(j)]
            else:
                g += d[A.index(j),A.index(A_i)]
        if g < g_min:
            g_min = g
            min_g_i = parti[rrr]            
    return min_g_i



###################################
def region_share_border_area(parti, area, neighbors):
    
    region_i = []
    nei = neighbors[area]
    
    for i in range(len(parti)):
        region = parti[i]
        for aaa in region:
            if aaa in nei:
                region_i.append(i)
    region_i = np.unique(region_i)
    return region_i
    
    
    

###################################
def p_np_gen(P,A):
    new_p = []
    tot = len(A)
    for li in P:
        sub_p = [1 if i+1 in li else 0 for i in range(tot)]
        new_p.append(sub_p)

    p_np = np.array([np.array(l) for l in new_p])
    return p_np



###################################
def neighbor_region(p_np, neighbors, area, region): #find the neignbor regions of an area
    neighhbors_region = []
    neighhbors_area = neighbors[area]
    for i in neighhbors_area:
        if i in region:
            continue
        else:
            rrr = np.where(p_np[:,i-1] == 1)[0] ### generate p_np!!!
            neighhbors_region.append(rrr)
    neighhbors_region = list(np.unique(np.array(neighhbors_region)))
    
    return neighhbors_region



###################################
def region_area_adj_matrix(p_np, P, neighbors,A): #find the neignbor regions of a region
    region_area_adj = np.zeros((p_np.shape))
    for region in P:
        for area in region:
            neighhbors_region_ = neighbor_region(p_np, neighbors, area, region)
            for ii in range(len(neighhbors_region_)):
                region_area_adj[neighhbors_region_[ii],area-1] = 1
    return region_area_adj



###################################
def H(P,A,d):
    H = 0
    for region in P:
        h = 0
        for area_index_1 in range(len(region)): 
            for area_index_2 in range(area_index_1+1,len(region)):
                if A.index(region[area_index_1])<A.index(region[area_index_2]):
                    h += d[A.index(region[area_index_1]), A.index(region[area_index_2])]                   
                else:
                    h += d[A.index(region[area_index_2]), A.index(region[area_index_1])]
        H += h
    return H 



###################################
def H_np(p_reg_ind,d):
    H = sum(sum(np.multiply(p_reg_ind,d)))
    return H



###################################
def p_reg_ind_gen(p,A):
    p_reg_ind = np.zeros((len(A),len(A)))
    for region in p:
        for area_ind_1 in range(len(region)):
            for area_ind_2 in range(area_ind_1+1,len(region)): 
                if region[area_ind_1]<region[area_ind_2]:
                    p_reg_ind[region[area_ind_1]-1][region[area_ind_2]-1] = 1
                else:
                    p_reg_ind[region[area_ind_2]-1][region[area_ind_1]-1] = 1
    return p_reg_ind



###################################
def H_per_region(p,A,d):
    
    H_region = []
    for region in p:
        p_reg_ind_r = np.zeros((len(A),len(A)))
        for area_ind_1 in range(len(region)):
            for area_ind_2 in range(area_ind_1+1,len(region)): 
                if region[area_ind_1]<region[area_ind_2]:
                    p_reg_ind_r[region[area_ind_1]-1][region[area_ind_2]-1] = 1
                else:
                    p_reg_ind_r[region[area_ind_2]-1][region[area_ind_1]-1] = 1
        H_r = sum(sum(np.multiply(p_reg_ind_r,d)))
        H_region.append(H_r)

    return H_region



################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################



###################################
def MOE_perc_region_VRT_multi(region,A,data,data_vrt,count_data,z):
    
    m = dict()
    MOE = dict()
    MOE_perc = dict()

    for iii in data:
        data_mean_ = data[iii]
        data_vrt_ = data_vrt[iii]
        count_data_ = count_data[iii]

        vr_data_sum = np.zeros((1,len(data_vrt_[0])))
        m_ = 0

        for i in region:
            ind_i = A.index(i)
            vr_data_sum = vr_data_sum + data_vrt_[ind_i]
            m_ += data_mean_[ind_i]
        if count_data_ == 0:
            m_ = m_/len(region)
        MOE_ = z*np.std(vr_data_sum)
        MOE_perc_ = MOE_/m_

        m[iii] = m_
        MOE[iii] = MOE_
        MOE_perc[iii] = MOE_perc_
              
    return m, MOE, MOE_perc



###################################
def GrowRegions_vr_multi(A,neighbors,data,d,data_vrt,data_moe,threshold,tighter_param,z,count_data, PRINT):
    
    tighter_threshold = threshold * tighter_param
    if PRINT == 1:
        print("Tighter threshold:",threshold,"*",tighter_param,"=",tighter_threshold) 

    # Grow regions from initial seeds
    # such that the value of attribute l in each region is above threshold.

    for_enclave_partitions = []
    e = []
    A_u = copy.deepcopy(A)
    A_a = []

    while A_u != []:
        k = rd.randrange(0, len(A_u), 1)
        A_k = A_u[k]
        A_u.pop(k)
        A_a.append(A_k)      
        ind_k = A.index(A_k)
        
        m = dict()         ########
        MOE = dict()       ########
        MOE_perc = dict()  ########        
        for iii in data:
            m[iii] = data[iii][ind_k]
            MOE[iii] = data_moe[iii][ind_k]
            MOE_perc[iii] = MOE[iii]/m[iii]

        threshold_hold = 1
        for iii in MOE_perc:
            if MOE_perc[iii] > tighter_threshold:
                threshold_hold = 0

        if threshold_hold == 1:
            R_k = [A_k]
            for_enclave_partitions.append(R_k)

        if threshold_hold == 0:
            R_k = [A_k]
            if PRINT == 1:
                print("Seed area: ", R_k)

            N = copy.deepcopy(neighbors[A_k])
            for aaa in A_a:
                if aaa in N:
                    N.pop(N.index(aaa))

            feasible = 1

            while threshold_hold == 0:
                if N!= []:
                    i = min_g_A(N, R_k, d, A)
                    if PRINT == 1:
                        print("Zone added: ",i)
                    R_k.append(i)
                    N.pop(N.index(i))
                    for aaa in neighbors[i]:
                        N.append(aaa)
                    for aaa in A_a:
                        if aaa in N:
                            N.pop(N.index(aaa))
                    N = list(np.unique(np.array(N)))
                    m, MOE, MOE_perc = MOE_perc_region_VRT_multi(R_k,A,data,data_vrt,count_data,z)

                    if PRINT == 1:
                        for iii in MOE_perc:
                            print("Data"+str(iii)+"MOE% = ",MOE_perc[iii]*100,"%")

                    if i not in A_u:
                        print('i not in A_u!!!')
                        print("i:",i)
                        print('A_u:',A_u)

                    A_u.pop(A_u.index(i))
                    A_a.append(i)

                threshold_hold = 1
                for iii in MOE_perc:
                    if MOE_perc[iii] > tighter_threshold:
                        threshold_hold = 0
                        if PRINT == 1:
                            print("Data"+str(iii)+"MOE%: NOT GOOD!")
                        break

                if N == [] and threshold_hold == 0:
                    for aaa in R_k:
                        e.append(aaa)

                    feasible = 0
                    if PRINT == 1:
                        print("Infeasible!!!")

                    for aaa in R_k:
                        A_a.pop(A_a.index(aaa))

                    if PRINT == 1:    
                        print("Unassigned:", A_u)
                    break

            if feasible == 1:
                R_k = list(np.sort(R_k))
                for_enclave_partitions.append(R_k)

            if PRINT == 1:
                print("for_enclave_partitions:",for_enclave_partitions)
                print("----------")
                
    return for_enclave_partitions, e, A_a



###################################
def GrowRegions_enclaveassign_example_vr_multi(maxitr, A,neighbors,data,d,data_vrt,data_moe,threshold,tighter_param, count_data,z,PRINT):

    feasi_partitions = []
    for_enclave_partitions = []
    maxP = 0

    #construction phase
    for i in range(1,maxitr):
        if PRINT == 1:
            print("--------------")
            print("iteration ",i)

        parti, e, A_a = GrowRegions_vr_multi(A,neighbors,data,d,data_vrt,data_moe,threshold,tighter_param,z, count_data, PRINT)

        if PRINT == 1:
            print("GrowRegion: ")
            print("Partition: ",parti)
            print("Enclave areas: ",e)
            print("Assigned areas: ", A_a)

        p = len(parti) #number of resgions in "parti"
        if p > maxP:
            for_enclave_partitions = [parti]
            maxP = p
        elif p == maxP:
            for_enclave_partitions.append(parti)
        elif p < maxP:
            pass

    for_enclave_partitions_uni = []
    for i in for_enclave_partitions:
        if i not in for_enclave_partitions_uni:
            for_enclave_partitions_uni.append(i)
    
    if PRINT == 1:
        print("-------------------------------")
        print("FINAL for Enclave Partitions: ")
        for i in range(len(for_enclave_partitions_uni)):
            print(for_enclave_partitions_uni[i])

    for parti in for_enclave_partitions_uni:
        P_feasi, e = AssignEnclaves(parti,d,A,neighbors) ######## func
        feasi_partitions.append(P_feasi)
        if PRINT == 1:
            print(P_feasi)
            print("Enclaves:", e)

    feasi_partitions_uni = []
    for i in feasi_partitions:
        if i not in feasi_partitions_uni:
            feasi_partitions_uni.append(i)

    if PRINT == 1:
        print("-------------------------------")
        print("P_feasi after enclave assign: ")
    for i in range(len(feasi_partitions)): 
        for ind_reg in range(len(feasi_partitions[i])):
            feasi_partitions[i][ind_reg] = list(np.sort(np.array(feasi_partitions[i][ind_reg])))
        if PRINT == 1:
            print(feasi_partitions[i])

    if PRINT == 1:
        print("------------")

    return feasi_partitions, for_enclave_partitions_uni, e



###################################
def check_feasibility_vr_multi(P,A,data_vrt,data,z,count_data,threshold):   
    for region in P:
        m, MOE, MOE_perc = MOE_perc_region_VRT_multi(region,A,data,data_vrt,count_data,z)

        threshold_hold = 1
        for iii in MOE_perc:
            if MOE_perc[iii] > threshold:
                threshold_hold = 0  
    return threshold_hold



###################################
def neighbor_set_1swap_vr_multi(p,A,data_vrt,data,d,z,threshold,count_data,neighbors,PRINT): 

    neighbor_sol = []
    neighbor_H = []

    ### generate the region_area_adj_matrix
    p_np = p_np_gen(p,A)
    r_a_adj = region_area_adj_matrix(p_np, p, neighbors,A)

    p_reg_ind_ori = p_reg_ind_gen(p,A)
    H_p = H_np(p_reg_ind_ori, d)
    #     print(p_reg_ind_ori)

    H_region = H_per_region(p,A,d)

    counttt = 0
    switch_list = []
    pair_weight_list = []

    for reg_ind_1 in range(len(p)):
        for reg_ind_2 in range(reg_ind_1+1,len(p)):
            switch_list.append((reg_ind_1,reg_ind_2))
            pair_weight = H_region[reg_ind_1] + H_region[reg_ind_2]
            pair_weight_list.append(pair_weight)

    new = [(i,j) for i,j in zip(switch_list, pair_weight_list)]
    sorted_ = sorted(new, key = lambda x:x[1])
    switch_list_sorted = [i[0] for i in sorted_]
    pair_weight_list_sorted = [i[1] for i in sorted_]

    for itr in range(len(switch_list_sorted)):
        pair_H_accu_list = []
        pair_H_accu = 0
        for iii in range(len(switch_list_sorted)):
            pair_H_accu += pair_weight_list_sorted[iii]
            pair_H_accu_list.append(pair_H_accu)
                
        k = rd.randrange(0, pair_H_accu_list[-1]+1, 1)
        pair = switch_list_sorted[np.where((np.array(pair_H_accu_list) >= k) == True)[0][0]]

        pop_ind = switch_list_sorted.index(pair)
        switch_list_sorted.pop(pop_ind)
        pair_weight_list_sorted.pop(pop_ind)

        reg_ind_1 = pair[0]
        reg_ind_2 = pair[1]

        aaa = p_np[reg_ind_1,:]
        bbb = r_a_adj[reg_ind_2,:]
        ind_a = np.where(aaa == 1)[0]
        ind_b = np.where(bbb == 1)[0]
        border_areas_1_2 = intersection(ind_a,ind_b)

        for aaa in border_areas_1_2:
            region_1_copy = copy.deepcopy(p[reg_ind_1])
            region_2_copy = copy.deepcopy(p[reg_ind_2])

            region_1_copy.pop(region_1_copy.index(aaa+1))
            region_2_copy.append(aaa+1)

            ### 连续性check
            reach_1 = 1
            for iii in range(1,len(region_1_copy)):
                if check_connection(region_1_copy[0],region_1_copy[iii],region_1_copy,neighbors) == 1:
                    pass
                else:
                    reach_1 = 0
                    if PRINT == 1:
                        print('region not connected: ',region_1_copy)
                    break
            reach_2 = 1
            for iii in range(1,len(region_2_copy)):
                if check_connection(region_2_copy[0],region_2_copy[iii],region_2_copy,neighbors) == 1:
                    pass
                else:
                    reach_2 = 0
                    if PRINT == 1:
                        print('region not connected: ',region_2_copy)
                    break

            if  reach_1 == 1 and reach_2 == 1: 
                new_nei = []
                for reg in p:
                    if reg == p[reg_ind_1] or reg == p[reg_ind_2]:
                        pass
                    else:
                        new_nei.append(reg)
                new_nei.append(region_1_copy)
                new_nei.append(region_2_copy)
    #                     print(new_nei)
                ### feasibility check
                feasi = check_feasibility_vr_multi(new_nei,A,data_vrt,data,z,count_data,threshold)
                if feasi == 1:
                    neighbor_sol.append(new_nei)
                    ## p_reg_ind_ori --> new_nei_reg_ind
                    new_nei_reg_ind = copy.deepcopy(p_reg_ind_ori)
                    for area in p[reg_ind_1]:
                        if area<(aaa+1):
                            new_nei_reg_ind[area-1][aaa] = 0
                        else:
                            new_nei_reg_ind[aaa][area-1] = 0
                    for area in p[reg_ind_2]:
                        if area<(aaa+1):
                            new_nei_reg_ind[area-1][aaa] = 1
                        else:
                            new_nei_reg_ind[aaa][area-1] = 1
                    neighbor_H.append(H_np(new_nei_reg_ind, d))
                    if H_np(new_nei_reg_ind, d) < H_p:
    #                             print("stopped",H_np(new_nei_reg_ind, d),H_p)
                        return [new_nei],[H_np(new_nei_reg_ind, d)]
    #                         print(new_nei_reg_ind)
    #                         print(H_np(new_nei_reg_ind, d))
                    if PRINT == 1:
                        print('Feasible!! Appended.')

        aaa = p_np[reg_ind_2,:]
        bbb = r_a_adj[reg_ind_1,:]
        ind_a = np.where(aaa == 1)[0]
        ind_b = np.where(bbb == 1)[0]
        border_areas_2_1 = intersection(ind_a,ind_b)

        for aaa in border_areas_2_1:
            region_2_copy = copy.deepcopy(p[reg_ind_2])
            region_1_copy = copy.deepcopy(p[reg_ind_1])

            region_2_copy.pop(region_2_copy.index(aaa+1))
            region_1_copy.append(aaa+1)

            ### 连续性check
            reach_1 = 1
            for iii in range(1,len(region_1_copy)):
                if check_connection(region_1_copy[0],region_1_copy[iii],region_1_copy,neighbors) == 1:
                    pass
                else:
                    reach_1 = 0
                    if PRINT == 1:
                        print('region not connected: ',region_1_copy)
                    break
            reach_2 = 1
            for iii in range(1,len(region_2_copy)):
                if check_connection(region_2_copy[0],region_2_copy[iii],region_2_copy,neighbors) == 1:
                    pass
                else:
                    reach_2 = 0
                    if PRINT == 1:
                        print('region not connected: ',region_2_copy)
                    break

            if  reach_1 == 1 and reach_2 == 1: 
                new_nei = []
                for reg in p:
                    if reg == p[reg_ind_1] or reg == p[reg_ind_2]:
                        pass
                    else:
                        new_nei.append(reg)
                new_nei.append(region_1_copy)
                new_nei.append(region_2_copy)
    #                     print(new_nei)
                ### feasibility check
                feasi = check_feasibility_vr_multi(new_nei,A,data_vrt,data,z,count_data,threshold)
                if feasi == 1:
                    neighbor_sol.append(new_nei)
                    ## p_reg_ind_ori --> new_nei_reg_ind
                    new_nei_reg_ind = copy.deepcopy(p_reg_ind_ori)
                    for area in p[reg_ind_2]:
                        if area<(aaa+1):
                            new_nei_reg_ind[area-1][aaa] = 0
                        else:
                            new_nei_reg_ind[aaa][area-1] = 0
                    for area in p[reg_ind_1]:
                        if area<(aaa+1):
                            new_nei_reg_ind[area-1][aaa] = 1
                        else:
                            new_nei_reg_ind[aaa][area-1] = 1
                    neighbor_H.append(H_np(new_nei_reg_ind, d))
                    if H_np(new_nei_reg_ind, d) < H_p:
    #                             print("stopped",H_np(new_nei_reg_ind, d),H_p)
                        return [new_nei],[H_np(new_nei_reg_ind, d)]
    #                         print(new_nei_reg_ind)
    #                         print(H_np(new_nei_reg_ind, d))
                    if PRINT == 1:
                        print('Feasible!! Appended.')

    ### rank them based on H (small -> large)
    neighbor_new = [(i,j) for i,j in zip(neighbor_sol, neighbor_H)]
    neighbor_sorted = sorted(neighbor_new, key = lambda x:x[1])
    neighbor_sol_sorted = [i[0] for i in neighbor_sorted]
    neighbor_H_sorted = [i[1] for i in neighbor_sorted]
    for row_sol_ind in range(len(neighbor_sol_sorted)):
        row_sol = neighbor_sol_sorted[row_sol_ind]
        row_avg = [np.mean(i) for i in row_sol]
        row_sol_new = [(i,j) for i,j in zip(row_sol, row_avg)]
        row_sorted = sorted(row_sol_new, key=lambda x:x[1])
        row_sol_sorted = [i[0] for i in row_sorted]
        neighbor_sol_sorted[row_sol_ind] = row_sol_sorted
    if PRINT == 1:
        print("Sorted neighbors and H:")
        for i in range(len(neighbor_sol_sorted)):
    #             print(neighbor_sol_sorted[i])
            print(neighbor_H_sorted[i])

    #     print("H computing of the neighbors finished!")


    return neighbor_sol_sorted, neighbor_H_sorted



###################################
def Tabu_search_vr_multi(Tabu_len, max_iter, feasi_partitions, A, data, d, data_vrt, neighbors, SWAP, z, threshold, count_data, PRINT):

    p_best_overall = feasi_partitions[0]
    
    plot_iter = []
    plot_time = []
    plot_time_cumu = []
    plot_H = []
      
    try:
        for ppp in range(len(feasi_partitions)):
            p = feasi_partitions[ppp]

            print(" ")
            print("Start from the Partition ",ppp," from GrowRegion")
            print("$$$$$$$$$$$$$$$$$")

            best_p_record = []
            best_p_x_axis = []

            Tabu_list = []
            Tabu_iter = []

            p_current = p
            p_best = p
            p_best_reg_ind = p_reg_ind_gen(p_best,A)
            H_p_best = H_np(p_best_reg_ind,d)        

            if PRINT == 1:
                print(p_best, A, d)
                print("H_best: ",H_p_best)
            Tabu_list.append(p_current)
            Tabu_iter.append(1)

            for iiii in range(max_iter): 
                start_iiii = time.time()
                if SWAP == 1:
                    neighbor_sol_sorted, neighbor_H_sorted = neighbor_set_1swap_vr_multi(p_current,A,data_vrt,data,d,z,threshold,count_data,neighbors,PRINT)
                else:
                    neighbor_sol_sorted, neighbor_H_sorted = neighbor_set_2swap_vr(p_current,A,data_vrt,data,z,threshold,count_data,neighbors,PRINT)

                for nei_ind in range(len(neighbor_sol_sorted)):
                    nei_current = neighbor_sol_sorted[nei_ind] ## pick the best one from the neighborhood
                    if PRINT == 1:
                        print("H of neighbor: ", neighbor_H_sorted[nei_ind], "H of BEST: ", H_p_best, p_best)
                    ### move to the top 1
                    if nei_current not in Tabu_list:
                        if neighbor_H_sorted[nei_ind] < H_p_best:
                            # move to it
                            if PRINT == 1:
                                print("Move to the neighbor", nei_ind, ": ", nei_current, "(Better sol, not in Tabu)")
                            p_current = nei_current
                            p_best = nei_current
                            p_best_reg_ind = p_reg_ind_gen(p_best,A)
                            H_p_best = H_np(p_best_reg_ind,d)
                            # add it to Tabulist & set Tabuiter = tabu_len
                            Tabu_list.append(nei_current)
                            Tabu_iter.append(1)

                            ### check Tabu list iterations  
                            Tabu_iter = np.array(Tabu_iter)
                            pop_index = list(np.where(Tabu_iter == Tabu_len)[0])
                            Tabu_list = np.delete(Tabu_list, pop_index, axis=0)
                            Tabu_iter = np.delete(Tabu_iter, pop_index, axis=0)
                            Tabu_list = list(Tabu_list)
                            for i in range(len(Tabu_list)):
                                Tabu_list[i] = list(Tabu_list[i])
                            Tabu_iter = list(Tabu_iter)
                            for iii in range(len(Tabu_iter)):
                                Tabu_iter[iii] = Tabu_iter[iii]+1

                            break

                        else:
                            # move to it
                            if PRINT == 1:
                                print("Move to the neighbor", nei_ind, ": ", nei_current, "(Worsen allowed, not in Tabu)")
                            p_current = nei_current
                            # add it to Tabulist
                            Tabu_list.append(nei_current)
                            Tabu_iter.append(1)

                            ### check Tabu list iterations  
                            Tabu_iter = np.array(Tabu_iter)
                            pop_index = list(np.where(Tabu_iter == Tabu_len)[0])
                            Tabu_list = np.delete(Tabu_list, pop_index, axis=0)
                            Tabu_iter = np.delete(Tabu_iter, pop_index, axis=0)
                            Tabu_list = list(Tabu_list)
                            for i in range(len(Tabu_list)):
                                Tabu_list[i] = list(Tabu_list[i])
                            Tabu_iter = list(Tabu_iter)
                            for iii in range(len(Tabu_iter)):
                                Tabu_iter[iii] = Tabu_iter[iii]+1

                            break

                    else: # if neighbor in Tabu list
                        if neighbor_H_sorted[nei_ind] < H_p_best:
                            # move to it
                            if PRINT == 1:
                                print("Move to the neighbor", nei_ind, ": ", nei_current, "(Better sol, in Tabu)")
                            p_current = nei_current
                            p_best = nei_current
                            p_best_reg_ind = p_reg_ind_gen(p_best,A)
                            H_p_best = H_np(p_best_reg_ind,d)
                            # remove from Tabulist                        
                            pop_index = Tabu_list.index(nei_current)
                            Tabu_list = np.delete(Tabu_list, pop_index, axis=0)
                            Tabu_iter = np.delete(Tabu_iter, pop_index, axis=0)
                            Tabu_list = list(Tabu_list)
                            for i in range(len(Tabu_list)):
                                Tabu_list[i] = list(Tabu_list[i])
                            Tabu_iter = list(Tabu_iter)
                            # add it to Tabulist & set Tabuiter = tabu_len
                            Tabu_list.append(nei_current)
                            Tabu_iter.append(1)

                            ### check Tabu list iterations  
                            Tabu_iter = np.array(Tabu_iter)
                            pop_index = list(np.where(Tabu_iter == Tabu_len)[0])
                            Tabu_list = np.delete(Tabu_list, pop_index, axis=0)
                            Tabu_iter = np.delete(Tabu_iter, pop_index, axis=0)
                            Tabu_list = list(Tabu_list)
                            for i in range(len(Tabu_list)):
                                Tabu_list[i] = list(Tabu_list[i])
                            Tabu_iter = list(Tabu_iter)
                            for iii in range(len(Tabu_iter)):
                                Tabu_iter[iii] = Tabu_iter[iii]+1
                            Tabu_iter = list(Tabu_iter)

                            break
                        else:
                            # NOT move to it
                            if PRINT == 1:
                                print("NOT MOVING to the neighbor", nei_ind, ": ", nei_current, "(Worsen, in Tabu)")

                    ### check Tabu list iterations  
                    Tabu_iter = np.array(Tabu_iter)
                    pop_index = list(np.where(Tabu_iter == Tabu_len)[0])
                    Tabu_list = np.delete(Tabu_list, pop_index, axis=0)
                    Tabu_iter = np.delete(Tabu_iter, pop_index, axis=0)
                    Tabu_list = list(Tabu_list)
                    for i in range(len(Tabu_list)):
                        Tabu_list[i] = list(Tabu_list[i])
                    Tabu_iter = list(Tabu_iter)
                    for iii in range(len(Tabu_iter)):
                        Tabu_iter[iii] = Tabu_iter[iii]+1

                best_p_record.append(H_p_best)
                best_p_x_axis.append(iiii)
                
                end_iiii = time.time()
                time_iiii = end_iiii - start_iiii
                
                plot_iter.append(iiii)
                plot_H.append(H_p_best)
                plot_time.append(time_iiii)
                plot_time_cumu.append(sum(plot_time))
                if iiii%100 == 0 and iiii!=0:
                    
                    print("Iteration",iiii)
                    print("Current best: ", H_p_best)
                    
                    plt.plot(plot_iter, plot_time)
                    plt.title("Computation Time of Each Iteration")
                    plt.xlabel("Iteration")
                    plt.ylabel("Seconds")
                    plt.savefig("time_fig/" + str(iiii) + "iter_time.pdf")
                    plt.show()                    
                    
                    plt.plot(plot_iter, plot_H)
                    plt.title("Heterogeneity of Each Iteration")
                    plt.xlabel("Iteration")
                    plt.savefig("time_fig/" + str(iiii) + "iter_H.pdf")
                    plt.show()
                    
                    plt.plot(plot_time_cumu, plot_H)
                    plt.title("Heterogeneity Change with Computation Time")
                    plt.xlabel("Cumulative Computation Time (Seconds)")
                    plt.savefig("time_fig/" + str(iiii) + "time_H.pdf")
                    plt.show()
                   
            
            p_best_overall_reg_ind = p_reg_ind_gen(p_best_overall,A)
            H_p_best_overall = H_np(p_best_overall_reg_ind,d)

            if H_p_best < H_p_best_overall:
                p_best_overall = p_best        


        return p_best_overall

    
    except:
        print("Keyboard interrupted")
        print("Current best overall: ", H(p_best_overall, A, d))
        print("Current best of the starting point: ", H_p_best)
        if H_p_best < H(p_best_overall, A, d):
            p_best_overall = p_best
            
        plt.plot(plot_iter, plot_time)
        plt.title("Computation Time of Each Iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Seconds")
        plt.savefig("time_fig/" + str(iiii) + "iter_time.pdf")
        plt.show()                    

        plt.plot(plot_iter, plot_H)
        plt.title("Heterogeneity of Each Iteration")
        plt.xlabel("Iteration")
        plt.savefig("time_fig/" + str(iiii) + "iter_H.pdf")
        plt.show()

        plt.plot(plot_time_cumu, plot_H)
        plt.title("Heterogeneity Change with Computation Time")
        plt.xlabel("Cumulative Computation Time (Seconds)")
        plt.savefig("time_fig/" + str(iiii) + "time_H.pdf")
        plt.show()
        
        return p_best_overall
