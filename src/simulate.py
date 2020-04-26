from src.system import System


def simulate(system: System, time: float):
    # time in (fs)

    state_list = []

    while system.t < time:
        print(system.t / time)
        system.step()
        state_list.append(system.make_data())

    return state_list
