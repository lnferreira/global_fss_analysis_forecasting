# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

source("fire_season_analysis.R")

x = sample(x = 1:10, size = 200, replace = T)
x[c(25:34,55:65,120:180)] = 0
seasons = estimate_fire_season(x, num_seasons = 4, warn = F)
plot(x, t="l", xlab="time", ylab="fire counts")
for (i in 1:length(seasons))
    lines(seasons[[i]], x[seasons[[i]]], col=2, t='o')
