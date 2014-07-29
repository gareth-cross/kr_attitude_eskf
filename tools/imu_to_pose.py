#!/usr/bin/env python
"""Convert IMU orientations to poses for use with RViz.
    
    @param ~imu sensor_msgs.Imu topic to subscribe to.
    @note: Publishes the topic 'pose', which is a geometry_msgs.PoseStamped.
"""

import rospy
import message_filters

from geometry_msgs.msg import PoseStamped
from sensor_msgs.msg import Imu

class ImuToPose():
    msub=None
    mpub=None
    
    def __init__(self):
        self.mpub = rospy.Publisher('~pose', PoseStamped, queue_size=1)
        
        self.msub = message_filters.Subscriber('~imu', Imu)
        self.msub.registerCallback(self.callback)

    def callback(self, msg):
        if type(msg) is Imu:
            pose = PoseStamped()
            pose.header.stamp = msg.header.stamp
            pose.header.frame_id = "0"
            pose.pose.orientation = msg.orientation
            pose.pose.position.x = pose.pose.position.y = pose.pose.position.z = 0
            self.mpub.publish(pose)
        
def main():
    rospy.init_node('imu_to_pose')
    imuToPose = ImuToPose()
    rospy.spin()

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        import traceback
        traceback.print_exc()
