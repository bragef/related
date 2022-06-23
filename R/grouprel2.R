#--------------------------#
# Function - grouprel2
#--------------------------#
# This function calculates the average observed relatedness
# within groups, then generates "expected" values using
# simulations, and then compares the two. It is modified from
# the original function because it allows for inbreeding.
# It takes 4 arguments:
# 1. A data frame of genotype data
# 2. The name of the relatedness estimator to use
# 3. The names of the groups to consider in the analyses
# 4. The number of simulated data sets to generate ("iterations")
#------------------------------------------------------------#

grouprel2 <- function(genotypes, estimatorname, usedgroups, iterations) {
    
    #--- Identify which column holds appropriate relatedness value ---#
    if (estimatorname == "trioml") {
        estimator = 5
    }
    if (estimatorname == "wang") {
        estimator = 6
    }
    if (estimatorname == "lynchli") {
        estimator = 7
    }
    if (estimatorname == "lynchrd") {
        estimator = 8
    }
    if (estimatorname == "ritland") {
        estimator = 9
    }
    if (estimatorname == "quellergt") {
        estimator = 10
    }
    if (estimatorname == "dyadml") {
        estimator = 11
    }
    
    #--- Calculate relatedness ---#
    if (estimatorname == "trioml") {
        relatives <- coancestry(genotypes, trioml = 1, allow.inbreeding = TRUE)
    }
    if (estimatorname == "wang") {
        relatives <- coancestry(genotypes, wang = 1, allow.inbreeding = TRUE)
    }
    if (estimatorname == "lynchli") {
        relatives <- coancestry(genotypes, lynchli = 1, allow.inbreeding = TRUE)
    }
    if (estimatorname == "lynchrd") {
        relatives <- coancestry(genotypes, lynchrd = 1, allow.inbreeding = TRUE)
    }
    if (estimatorname == "ritland") {
        relatives <- coancestry(genotypes, ritland = 1, allow.inbreeding = TRUE)
    }
    if (estimatorname == "quellergt") {
        relatives <- coancestry(genotypes, quellergt = 1, allow.inbreeding = TRUE)
    }
    if (estimatorname == "dyadml") {
        relatives <- coancestry(genotypes, dyadml = 1, allow.inbreeding = TRUE)
    }
    
    #--- Separate file with just relatedness data ---$
    rels <- relatives$relatedness
    
    #---------------------------------------#
    # Find the unique values of group names #
    #---------------------------------------#
    if (usedgroups[1] == "all") {
        groupsall <- 1:length(rels[, 1])
        
        for (i in 1:length(rels[, 1])) {
            groupsall[i] <- substr(rels[i, 2], 1, 2)
        }
        
        groups <- unique(groupsall)
    } else {
        groups <- usedgroups
    }
    
    
    #--------------------------------------------#
    # Calculate Average Within Group Relatedness #
    #--------------------------------------------#
    
    #--- Create list of within-group comparison names ---#
    within <- paste(groups, groups, sep = "")
    
    #--- Find relevant relatedness values ---#
    #--- This loop takes a while!!        ---#
    relvalues <- 1:length(within)
    sizes <- 1:length(within)
    
    cat("\n Calculating within-group r-values...\n")
    for (i in 1:length(within)) {
        holder <- 0
        counter1 <- 0
        
        for (j in 1:length(rels[, 1])) {
            if (rels[j, 4] == within[i]) {
                holder <- holder + rels[j, estimator]
                counter1 <- counter1 + 1
            }
            
        }
        relvalues[i] <- holder / counter1
        sizes[i] <- counter1
        cat(sprintf("Group %s \t %f\n", within[i], relvalues[i]))
    }
    overallobs <- sum(relvalues) / length(relvalues)
    cat(sprintf("Overall \t %f\n", overallobs))
    
    #--- Write results to file ---#
    obsrel <- cbind(within, relvalues)
    write.csv(obsrel, "observed-r.csv")
    
    
    #---------------------------------------------#
    # Calculate Expected Within-Group Relatedness #
    #---------------------------------------------#
    
    #--- Create a matrix for holding all results ---#
    simresults <- data.frame(matrix(nrow = iterations, ncol = (length(within) + 1)))
    
    for (j in 1:iterations) {
        
        cat(sprintf("Iteration %d\n", j))
        
        #--- Create a randomized list for shuffling ---#
        randlist <- 1:length(genotypes[, 1])
        randlist <- sample(randlist, length(randlist), replace = FALSE)
        
        randgenos <- data.frame(matrix(nrow = length(randlist), ncol = length(genotypes[1, ])))
        
        #--- Re-organize genotype file ---#
        for (i in 1:length(randlist)) {
            randgenos[i, ] <- genotypes[randlist[i], ]
        }
        
        #--- Calculate Relatedness ---#
        if (estimatorname == "trioml") {
            simrels <- coancestry(randgenos, trioml = 1, allow.inbreeding = TRUE)
        }
        if (estimatorname == "wang") {
            simrels <- coancestry(randgenos, wang = 1, allow.inbreeding = TRUE)
        }
        if (estimatorname == "lynchli") {
            simrels <- coancestry(randgenos, lynchli = 1, allow.inbreeding = TRUE)
        }
        if (estimatorname == "lynchrd") {
            simrels <- coancestry(randgenos, lynchrd = 1, allow.inbreeding = TRUE)
        }
        if (estimatorname == "ritland") {
            simrels <- coancestry(randgenos, ritland = 1, allow.inbreeding = TRUE)
        }
        if (estimatorname == "quellergt") {
            simrels <- coancestry(randgenos, quellergt = 1, allow.inbreeding = TRUE)
        }
        if (estimatorname == "dyadml") {
            simrels <- coancestry(randgenos, dyadml = 1, allow.inbreeding = TRUE)
        }
        
        #--- Organize results into appropriately sized groups ---#
        counter1 <- 1
        counter2 <- 0
        
        for (k in 1:length(within)) {
            
            holder <- 0
            counter3 <- 0
            
            # Find size of group
            num <- sizes[k]
            
            # Calculate number of pairwise comparisons
            counter2 <- num
            
            for (l in counter1:(counter2 + counter1 - 1)) {
                holder <- holder + simrels$relatedness[l, estimator]
                counter3 <- counter3 + 1
            }
            simresults[j, k] <- holder / counter3
            counter1 <- counter2 + counter1
        }
        simresults[j, length(within) + 1] <- sum(simresults[j, 1:length(within)]) / length(within)
        
    }
    
    #--- Write results to file ---#
    write.csv(simresults, "expectedrel.csv")
    
    minx <- 0
    maxx <- 0
    #--- Plot results for each group ---#
    for (k in 1:length(within)) {
        
        if (min(simresults[, k]) < relvalues[k]) {
            minx <- min(simresults[, k]) - 0.2
        } else {
            minx <- relvalues[k] - 0.2
        }
        
        if (max(simresults[, k]) > relvalues[k]) {
            maxx <- max(simresults[, k]) + 0.2
        } else {
            maxx <- relvalues[k] + 0.2
        }
        
        
        hist(simresults[, k], main = within[k], xlim = c(minx, maxx), xlab = "Relatedness")
        arrows(x0 = relvalues[k], y0 = iterations * 0.15, x1 = relvalues[k], y1 = 0, col = "red", lwd = 3)
        ptest <- signif(((sum(simresults[, k] >= relvalues[k]) + 1) / iterations), 3)
        mtext(bquote(p < .(ptest)), side = 3)
    }
    
    #--- Plot overall results ---#
    if (min(simresults[ , length(within) + 1]) < overallobs) {
        minx <- min(simresults[ , length(within) + 1]) - 0.2
    } else {
        minx <- overallobs - 0.2
    }
    
    if (max(simresults[ , length(within) + 1]) > overallobs) {
        maxx <- max(simresults[ , length(within) + 1])  + 0.2
    } else {
        maxx <- overallobs + 0.2
    }
    
    hist(simresults[, length(within) + 1], main = "Overall", xlim = c(minx, maxx), xlab = "Relatedness")
    arrows(x0 = overallobs, y0 = iterations * 0.15, x1 = overallobs, y1 = 0, col = "red", lwd = 3)
    ptest <- signif(((sum(simresults[, length(within) + 1] >= overallobs) + 1) / iterations), 3)
    mtext(bquote(p < .(ptest)), side = 3)
}

