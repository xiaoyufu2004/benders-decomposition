import numpy as np
import gurobipy as gp
from gurobipy import GRB

# ------------------------------------------------------------------
# 数据：3 个仓库、4 个客户
# ------------------------------------------------------------------
F = np.array([100, 120, 130], dtype=float)     # 仓库固定成本
A = np.array([70, 60,  80],  dtype=float)      # 仓库产能
b = np.array([40, 30,  20, 20], dtype=float)   # 客户需求

d = np.array([[4, 6, 9, 7],                    # 运输成本 (3 × 4)
              [5, 4, 7, 6],
              [6, 3, 4, 5]], dtype=float)

warehouses = range(3)
customers  = range(4)
total_demand = b.sum()

# ------------------------------------------------------------------
# 主问题：x_i（二进制开仓变量） + θ（递归成本）
# ------------------------------------------------------------------
master = gp.Model("master")
master.Params.LogToConsole = 0

x = master.addVars(warehouses, vtype=GRB.BINARY, name="x")
theta = master.addVar(lb=0, name="theta")

# 目标：固定成本 + θ
master.setObjective(gp.quicksum(F[i] * x[i] for i in warehouses) + theta,
                    GRB.MINIMIZE)

# *可选* 约束：保证子问题可行（若想测试可行割，可注释掉此行）
# master.addConstr(gp.quicksum(A[i] * x[i] for i in warehouses) >= total_demand,
#                  name="basic_capacity")

# ------------------------------------------------------------------
# Benders 主循环
# ------------------------------------------------------------------
UB, LB = float("inf"), -float("inf")
EPS, MAX_ITER = 1e-6, 100

x_current = np.zeros(len(warehouses))            # 初始全开仓
best_x, best_y = None, None

for it in range(1, MAX_ITER + 1):
    print(f"\n=== 迭代 {it} ===")

    # ---------------- 子问题 (运输) ----------------
    sub = gp.Model("sub")
    sub.Params.LogToConsole = 0

    y = sub.addVars(warehouses, customers, lb=0, name="y")
    sub.setObjective(gp.quicksum(d[i, j] * y[i, j]
                                 for i in warehouses for j in customers),
                     GRB.MINIMIZE)

    # 产能 ≤ A_i * x̄_i
    for i in warehouses:
        sub.addConstr(gp.quicksum(y[i, j] for j in customers)
                      <= A[i] * x_current[i], name=f"cap_{i}")

    # 需求 = b_j
    for j in customers:
        sub.addConstr(gp.quicksum(y[i, j] for i in warehouses)
                      == b[j], name=f"dem_{j}")

    sub.optimize()

    # ---------- 1) 子问题可行 → 最优割 ----------
    if sub.Status == GRB.OPTIMAL:
        transp = sub.ObjVal
        total  = F @ x_current + transp

        if total < UB:
            UB, best_x = total, x_current.copy()
            best_y = np.array([[y[i, j].X for j in customers] for i in warehouses])

        print(f"子问题最优 | 运输 {transp:.2f} | 总 {total:.2f} | UB {UB:.2f}")

        cap_dual = np.array([sub.getConstrByName(f"cap_{i}").Pi for i in warehouses])
        dem_dual = np.array([sub.getConstrByName(f"dem_{j}").Pi for j in customers])

        rhs_const = dem_dual @ b                               # π_dem ᵀ b
        rhs_xpart = gp.quicksum(A[i] * cap_dual[i] * x[i]      # Σ A_i π_cap x_i
                                for i in warehouses)

        master.addConstr(theta >= rhs_const + rhs_xpart,
                         name=f"optcut_{it}")

    # ---------- 2) 子问题不可行 → 可行割 ----------
    elif sub.Status == GRB.INFEASIBLE:
        print("子问题不可行，生成可行割")
        sub.Params.InfUnbdInfo = 1
        sub.optimize()                                         # 触发 Farkas 信息

        cap_ray = np.array([sub.getConstrByName(f"cap_{i}").FarkasDual
                            for i in warehouses])
        dem_ray = np.array([sub.getConstrByName(f"dem_{j}").FarkasDual
                            for j in customers])

        expr = (dem_ray @ b) + gp.quicksum(A[i] * cap_ray[i] * x[i]
                                           for i in warehouses)

        # 计算当前 x̄ 的值以判定不等式方向
        viol = (dem_ray @ b) + np.sum(A * cap_ray * x_current)

        if viol > 0:
            master.addConstr(expr <= 0, name=f"feascut_{it}")  # πᵀ rhs > 0
            sign = "<="
        else:
            master.addConstr(expr >= 0, name=f"feascut_{it}")  # πᵀ rhs < 0
            sign = ">="
        print(f"添加可行割: πᵀrhs({viol:+.2f}) {sign} 0")

    else:
        raise RuntimeError(f"子问题状态异常：{sub.Status}")

    # ---------------- 求解主问题 ----------------
    master.optimize()
    if master.Status != GRB.OPTIMAL:
        raise RuntimeError("主问题求解失败")

    LB = master.ObjVal
    x_current = np.array([x[i].X for i in warehouses])
    gap = UB - LB
    print(f"主问题 | LB {LB:.2f} | UB {UB:.2f} | Gap {gap:.6f}")

    if gap <= EPS:
        print("\n*** 收敛 ***")
        break

# ------------------------------------------------------------------
# 输出结果
# ------------------------------------------------------------------
print("\n========== 最终结果 ==========")
print(f"最优目标值: {UB:.2f}\n")

print("仓库启用情况:")
for i in warehouses:
    print(f"  仓库 {i+1}: {'启用' if best_x[i] > 0.5 else '关闭'}")

print("\n运输矩阵 (仓库行 × 客户列):")
print(best_y)

print(f"\n固定成本 : {F @ best_x:.2f}")
print(f"运输成本 : {np.sum(d * best_y):.2f}")
print(f"总成本   : {UB:.2f}")
