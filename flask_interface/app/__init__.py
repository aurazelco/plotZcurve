# imports flask module
from flask import Flask

# defines app name
app=Flask(__name__)

# imports routes.py, which contains the python flask code to be run
from app import routes
