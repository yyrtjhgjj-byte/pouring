"""pourphys — ウォーターサーバーでコップに水を注ぐときの飛沫を見積もる半経験物理モデル。"""
from .props import Fluid, COLD, ROOM, HOT, PRESETS, G
from .server import Server
from .cup import Cup, Pose, place_cup
from .strategy import Strategy, Scenario, Evaluation, evaluate, evaluate_pour

__all__ = ["Fluid", "COLD", "ROOM", "HOT", "PRESETS", "G", "Server", "Cup", "Pose",
           "place_cup", "Strategy", "Scenario", "Evaluation", "evaluate", "evaluate_pour"]
