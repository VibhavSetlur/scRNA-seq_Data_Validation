#!/bin/bash

# Launch Script for snRNA-seq Pipeline Shiny App
# This script starts the Shiny app and opens it in your browser

set -e

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Launching snRNA-seq Pipeline Shiny App ===${NC}"
echo ""

# Check if app is already running
if curl -s http://localhost:3838 > /dev/null 2>&1; then
    echo -e "${YELLOW}⚠ Shiny app is already running!${NC}"
    echo -e "${GREEN}Opening in browser...${NC}"
    
    # Open in browser
    if command -v xdg-open &> /dev/null; then
        xdg-open http://localhost:3838
    elif command -v open &> /dev/null; then
        open http://localhost:3838
    else
        echo -e "${YELLOW}Please open your browser and go to: http://localhost:3838${NC}"
    fi
    
    echo ""
    echo -e "${GREEN}✅ App is ready at: http://localhost:3838${NC}"
    echo -e "${YELLOW}Press Ctrl+C to stop the app${NC}"
    exit 0
fi

# Start the Shiny app
echo -e "${BLUE}Starting Shiny app...${NC}"
echo -e "${YELLOW}This may take a few seconds...${NC}"
echo ""

# Start the app in background
Rscript scripts/run_shiny_app.R &
APP_PID=$!

# Wait for app to start
echo -e "${BLUE}Waiting for app to start...${NC}"
for i in {1..30}; do
    if curl -s http://localhost:3838 > /dev/null 2>&1; then
        echo -e "${GREEN}✅ App started successfully!${NC}"
        break
    fi
    sleep 1
    echo -n "."
done

echo ""
echo -e "${GREEN}Opening in browser...${NC}"

# Open in browser
if command -v xdg-open &> /dev/null; then
    xdg-open http://localhost:3838
elif command -v open &> /dev/null; then
    open http://localhost:3838
else
    echo -e "${YELLOW}Please open your browser and go to: http://localhost:3838${NC}"
fi

echo ""
echo -e "${GREEN}🎉 snRNA-seq Pipeline is ready!${NC}"
echo -e "${BLUE}Access the app at: http://localhost:3838${NC}"
echo ""
echo -e "${YELLOW}Press Ctrl+C to stop the app${NC}"

# Wait for user to stop
trap "echo ''; echo -e '${BLUE}Stopping Shiny app...${NC}'; kill $APP_PID 2>/dev/null; echo -e '${GREEN}App stopped.${NC}'; exit 0" INT

# Keep script running
wait $APP_PID
